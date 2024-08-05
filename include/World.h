//
//  This file is part of the FFEA simulation package
//
//  Copyright (c) by the Theory and Development FFEA teams,
//  as they appear in the README.md file.
//
//  FFEA is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  FFEA is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
//
//  To help us fund FFEA development, we humbly ask that you cite
//  the research papers on the package.
//

#ifndef WORLD_H_INCLUDED
#define WORLD_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <ctime>
#include <algorithm>

#include <boost/algorithm/string.hpp>
#include <typeinfo>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#ifdef FFEA_PARALLEL_FUTURE
#include <future>
#include <chrono>
#endif

// #include "MersenneTwister.h"
#include "RngStream.h"
#include "NearestNeighbourLinkedListCube.h"
#include "BEM_Poisson_Boltzmann.h"
#include "BiCGSTAB_solver.h"
#include "FFEA_user_info.h"
#include "FFEA_return_codes.h"
#include "FFEA_input_reader.h"
#include "mat_vec_types.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseSubstitutionSolver.h"
#include "Face.h"
#include "Blob.h"
#include "World.h"
#include "VdW_solver.h"
#include "Steric_solver.h"
#include "LJSteric_solver.h"
#include "GenSoftSSINT_solver.h"
#include "PreComp_solver.h"
#include "LJ_matrix.h"
#include "BindingSite.h"
#include "Spring.h"
#include "SparseMatrixFixedPattern.h"
#include "KineticState.h"
#include "rod_structure.h"
#include "rod_blob_interface.h"

#include "dimensions.h"
using namespace std;

class World {
friend struct ffea_test; //allow  our unit test class to see ffea_world's private
public:
    World();

    ~World();

    /* */
    int init(string FFEA_script_filename, int frames_to_delete, int mode, bool writeEnergy);

    /* */
    int get_smallest_time_constants();

    /* */
    int lem(set<int> blob_indices, int num_modes);

    /* */
    int dmm(set<int> blob_indices, int num_modes);

    /* */
    int dmm_rp(set<int> blob_indices, int num_modes);

    /* */
    int run();

    /* */
    int read_and_build_system(vector<string> script_vector);

    /* */
    int load_kinetic_maps(vector<string> map_fnames, vector<int> map_from, vector<int> map_to, int blob_index);

    /* */
    int build_kinetic_identity_maps();

    /* */
    int load_kinetic_states(string states_fname, int blob_index);

    /* */
    int load_kinetic_rates(string rates_fname, int blob_index);

    /* */
    void print_kinetic_rates_to_screen(int type);

    /* */
    void get_system_CoM(vector3 *system_CoM);

    /* */
    void get_system_centroid(vector3 *centroid);

    /* */
    void get_system_dimensions(vector3 *dimenstion_vector);

    /* */
    int enm(int *blob_index, int num_modes);

    /* */
    int get_num_blobs();


private:

    /** @brief 2-D Array of Blob objects (blob i, conformation j) */
    Blob **blob_array;

    /** @brief Which conformation is active in each blob */
    Blob **active_blob_array;

    /** @brief 1-D array containing pointers to all rod objects */
    rod::Rod **rod_array;

    /** @brief 1-D array containing pointers to all rod-blob interfaces */
    rod::Rod_blob_interface **rod_blob_interface_array;

    /** @brief Maps for kinetic switching of conformations */
    SparseMatrixFixedPattern ***kinetic_map;
    SparseMatrixFixedPattern ****kinetic_return_map;

    //@{
    /** @brief Kinetic State and Rate objects */
    KineticState **kinetic_state;
    scalar ***kinetic_rate;
    scalar ***kinetic_base_rate;
    //@}

    /** @brief An array of springs which connect nodes if necessary */
    Spring *spring_array;

    /** @brief And how many springs are there? */
    int num_springs;

    /** @brief How many kinetic binding sites are there? */
    int num_binding_sites;

    /** @brief Check whether mass is present anywhere, to determine whether or not to write kinetic energies to files */
    bool mass_in_system;

    /** @brief How many threads are available for parallelisation */
    int num_threads;

    /** @brief An array of pointers to random number generators (for use in parallel) */
    RngStream *rng;

    /** @brief A pointer to an array of arrays, containing the seeds of the different RNGStreams */
    unsigned long **Seeds;

    /** @brief The number of seeds stored in Seeds. */
    int num_seeds;

    /** @brief An array of pointers to random number generators for use in kinetics */
    RngStream *kinetic_rng;

    /** @brief Parameters being used for this simulation */
    SimulationParams params;

    /** @brief
     * Data structure keeping track of which `cell' each face lies in (where the world has been discretised into a grid of cells of dimension 1.5 kappa)
     * so that the BEM matrices may be constructed quickly and sparsely.
     */
    NearestNeighbourLinkedListCube lookup;

    /** @brief Output trajectory file */
    FILE *trajectory_out;

    /** @brief Output kinetics file */
    FILE *kinetics_out;

    /** @brief Output measurement file */
    FILE *measurement_out;

    /** @brief Output detailed measurements file. May be unneccesary */
    bool writeDetailed;
    FILE *detailed_meas_out;

    /** @brief Output file for the trajectory beads. Completely optional. */
    FILE *trajbeads_out;

    /** Reader objects */
    FFEA_input_reader *ffeareader;
    FFEA_input_reader *systemreader;

    //@{
    /** Energies */
    scalar kineticenergy, strainenergy, springenergy, **springfieldenergy, ssintenergy, preCompenergy;
    //@}

    /** Momenta */
    vector3 L;

    //@{
    /** Geometries */
    vector3 CoM, CoG;
    scalar rmsd;
    //@}

    /** @brief * Output Checkpoint file */
    FILE *checkpoint_out;

    /*
     *
     */
    //		SurfaceElementLookup surface_element_lookup;

    /** @brief BEM solver for the exterior electrostatics */
    BEM_Poisson_Boltzmann PB_solver;

    /** @brief Number of surface faces in entire system */
    int total_num_surface_faces;

    /** @brief
     * Vector of the electrostatic potential on each surface in entire system
     */
    scalar *phi_Gamma;
    scalar *J_Gamma;
    scalar *work_vec;

    /** @brief
     * Biconjugate gradient stabilised solver for nonsymmetric matrices
     */
    BiCGSTAB_solver nonsymmetric_solver;

    /** Van der Waals solver */
    VdW_solver *vdw_solver;

    /** @brief LJ parameters matrix */
    SSINT_matrix ssint_matrix;
    SSINT_matrix rod_lj_matrix;

    /** @brief Binding Interactions matrix */
    BindingSite_matrix binding_matrix;


    /**
      * @brief stores info within the <precomp> block at the .ffea file.
      */
    PreComp_params pc_params;
    /**
      * @brief PreComputed potentials solver
      */
    PreComp_solver pc_solver;


    vector3 box_dim;

    long long step_initial;

    int load_springs(const char *fname);

    rod::Rod_blob_interface* rod_blob_interface_from_block(vector<string> block, int interface_id, FFEA_input_reader* systemreader, rod::Rod** rod_array, Blob** blob_array);

    rod::Rod* rod_from_block(vector<string> block, int block_id, FFEA_input_reader* systemreader);

    void update_rod_steric_nbr_lists(rod::Rod* rod_a, rod::Rod* rod_b);

    void update_rod_vdw_nbr_lists(rod::Rod *rod_a, rod::Rod *rod_b, SSINT_matrix *lj_matrix);

    void rod_pbc_wrap(rod::Rod* current_rod, std::vector<float> dim);

    void rod_box_length_check(rod::Rod *current_rod, std::vector<float> dim);

    void activate_springs();

    int apply_springs();

    scalar get_spring_field_energy(int index0, int index1);

    /** @brief calculates the kinetic rates as a function of the energy of the system*/
    int calculate_kinetic_rates();

    /** @brief randomly chooses a new kinetic state based upon the kinetic rates / switching probabilitie */
    int choose_new_kinetic_state(int blob_index, int *target);

    /** @brief changes the kinetic state based upon the kinetic rates. Maps between conformations and adds/ removes bound sites */
    int change_kinetic_state(int blob_index, int target_state);

    int get_next_script_tag(FILE *in, char *buf);

    void apply_dense_matrix(scalar *y, scalar *M, scalar *x, int N);

    void do_es();

    void make_trajectory_from_eigenvector(string traj_out_fname, int blob_index, int mode_index, Eigen_VectorX evec, scalar step);

    void print_evecs_to_file(string fname, Eigen_MatrixX ev, int num_rows, int num_modes);

    void print_evals_to_file(string fname, Eigen_VectorX ev, int num_modes, scalar scale);

    void write_eig_to_files(scalar *evals_ordered, scalar **evecs_ordered, int num_modes, int num_nodes);

    void write_output_header(FILE *fout, string fname);

    void print_trajectory_and_measurement_files(int step, scalar wtime);
    void print_checkpoints();
    void write_pre_print_to_trajfile(int step);
    void do_nothing();

    int prebuild_nearest_neighbour_lookup_wrapper(scalar cell_size);
#ifdef FFEA_PARALLEL_FUTURE
    std::future<void> thread_writingTraj;
    std::future<int> thread_updatingVdWLL;
    std::future<int> thread_updatingPCLL;
    bool updatingVdWLL(); ///< check if the thread has been catched.
    bool updatingVdWLL_ready_to_swap(); ///< true if thread waiting to be catched.
    int catch_thread_updatingVdWLL(int step, scalar wtime, int where);
    bool updatingPCLL(); ///< check if the thread has been catched.
    bool updatingPCLL_ready_to_swap(); ///< true if thread waiting to be catched.
    int catch_thread_updatingPCLL(int step, scalar wtime, int where);
#endif

    void make_measurements();

    void write_measurements_to_file(FILE *fout, int step);

    void write_detailed_measurements_to_file(FILE *fout);

    void print_trajectory_conformation_changes(FILE *fout, int step, int *from_index, int *to_index);

    void print_kinetic_files(int step);

    void print_static_trajectory(int step, scalar wtime, int blob_index);

    /** @brief calculates the blob to blob corrections due to periodic boundary conditions*/
    void calc_blob_corr_matrix(int num_blobs,scalar *blob_corr);

    scalar *blob_corr;

    int die_with_dignity(int step, scalar wtime);
};

#endif
