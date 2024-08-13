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

#include "World.h"
#ifdef USE_MPI
#include "mpi.h"
#endif
#include <filesystem>


World::World()
{

    // Initialise everything to zero
    blob_array = NULL;
    rod_array = NULL;
    spring_array = NULL;
    kinetic_map = NULL;
    kinetic_return_map = NULL;
    kinetic_state = NULL;
    kinetic_rate = NULL;
    kinetic_base_rate = NULL;
    num_springs = 0;
    mass_in_system = false;
    num_threads = 1;
    rng = NULL;
    phi_Gamma = NULL;
    total_num_surface_faces = 0;
    box_dim.x = 0;
    box_dim.y = 0;
    box_dim.z = 0;
    step_initial = 0;
    trajectory_out = NULL;
    measurement_out = NULL;
    detailed_meas_out = NULL;
    writeDetailed = true;
    kinetics_out = NULL;
    checkpoint_out = NULL;
    vdw_solver = NULL;
    Seeds = NULL;
    num_seeds = 0;
    blob_corr = NULL;

    ffeareader = NULL;
    systemreader = NULL;
    kineticenergy = 0.0;
    strainenergy = 0.0;
    springenergy = 0.0;
    springfieldenergy = NULL;
    ssintenergy = 0.0;
    preCompenergy = 0.0;

    vector3_set_zero(L);
    vector3_set_zero(CoM);
    vector3_set_zero(CoG);
    rmsd = 0.0;
}

World::~World()
{
    delete[] rng;
    rng = NULL;
    num_threads = 0;

    for (int i = 0; i < params.num_blobs; ++i)
    {
        delete[] blob_array[i];
    }
    delete[] blob_array;
    blob_array = NULL;
    delete[] active_blob_array;
    active_blob_array = NULL;
    delete[] rod_array;
    rod_array = NULL;

    delete[] spring_array;
    spring_array = NULL;
    num_springs = 0;

    mass_in_system = false;

    delete[] kinetic_map;
    kinetic_map = NULL;

    delete[] kinetic_return_map;
    kinetic_return_map = NULL;

    delete[] kinetic_state;
    kinetic_state = NULL;
    delete[] kinetic_rate;
    kinetic_rate = NULL;
    delete[] kinetic_base_rate;
    kinetic_base_rate = NULL;

    delete[] phi_Gamma;
    phi_Gamma = NULL;

    delete[] blob_corr;

    for (int i = 0; i < num_seeds; i++)
    {
        delete[] Seeds[i];
    }
    delete[] Seeds;
    Seeds = NULL;
    num_seeds = 0;

    total_num_surface_faces = 0;

    box_dim.x = 0;
    box_dim.y = 0;
    box_dim.z = 0;
    step_initial = 0;

    if (trajectory_out != NULL)
    {
        fclose(trajectory_out);
    }
    if (trajectory_out != NULL)
    {
        fclose(measurement_out);
    }
    if (trajectory_out != NULL)
    {
        fclose(checkpoint_out);
    }
    if ((writeDetailed) && (detailed_meas_out != NULL))
    {
        fclose(detailed_meas_out);
    }
    writeDetailed = false;
    if (params.calc_kinetics == 1)
    {
        fclose(kinetics_out);
    }
    trajectory_out = NULL;
    measurement_out = NULL;
    checkpoint_out = NULL;
    detailed_meas_out = NULL;
    kinetics_out = NULL;

    delete vdw_solver;
    vdw_solver = NULL;

    delete ffeareader;
    delete systemreader;
    ffeareader = NULL;
    systemreader = NULL;
    kineticenergy = 0.0;
    strainenergy = 0.0;
    springenergy = 0.0;

    delete[] springfieldenergy;
    springfieldenergy = NULL;

    ssintenergy = 0.0;
    preCompenergy = 0.0;

    vector3_set_zero(L);
    vector3_set_zero(CoM);
    vector3_set_zero(CoG);
    rmsd = 0.0;
}

/**
 * @brief Reads the .ffea file and initialises the World.
 * @param[in] string FFEA_script_filename
 * @details Open and read .ffea file,
 *   parse the <param> block through SimulationParams::extract_params in the
 *   "SimulationParams params" private attribute,
 *   parse the <blobs> and <springs> blocks through World::read_and_build_system
 * initialise a number of RNG,
 * prepare output files,
 * initialise VdW solver,
 * initialise BEM PBE solver
 * */

int World::init(string FFEA_script_filename, int frames_to_delete, int mode, bool writeDetail)
{

    // Set some constants and variables
    this->writeDetailed = writeDetail;
    int i, j, k;

    string buf_string;
    ffeareader = new FFEA_input_reader();
#ifdef USE_MPI
    double st, et;

    st = MPI::Wtime();
#endif

    // Copy entire script into string
    vector<string> script_vector;
    if (ffeareader->file_to_lines(FFEA_script_filename, &script_vector) == FFEA_ERROR)
    {
        return FFEA_ERROR;
    }

    // Get params section
    cout << "Extracting Parameters..." << endl;
    params.FFEA_script_filename = FFEA_script_filename; // includes absolute path.
    if (params.extract_params(script_vector) != 0)
    {
        FFEA_error_text();
        printf("Error parsing parameters in SimulationParams::extract_params()\n");
        return FFEA_ERROR;
    }
    cout << "...done!" << endl;

    // Check for consistency
    cout << "\nVerifying Parameters..." << endl;
    if (params.validate(mode) != 0)
    {
        FFEA_error_text();
        printf("Parameters found to be inconsistent in SimulationParams::validate()\n");
        return FFEA_ERROR;
    }
    if ((params.num_blobs) == 1)
    {
        writeDetailed = false;
        printf("\n\tA single blob is simulated, and thus the detailed measurements would be redundant and are not needed\n");
    }

    // Build kinetic maps if necessary (and rates and binding site matrix). These are at the World level in case global kinetic calculations are ever included
    if (params.calc_kinetics == 1)
    {

        // Load the binding params matrix
        if (binding_matrix.init(params.bsite_in_fname) == FFEA_ERROR)
        {
            FFEA_ERROR_MESSG("Error when reading from binding site params file.\n")
        }

        // A 3D matrix describing the switching rates for each blob i.e. kinetic_rate[blob_index][from_conf][to_conf]
        kinetic_rate = new scalar **[params.num_blobs];
        kinetic_base_rate = new scalar **[params.num_blobs];

        // A 3D matrix holding the position maps enabline the switch from one conformation to another i.e. kinetic_map[blob_index][from_conf][to_conf]
        kinetic_map = new SparseMatrixFixedPattern **[params.num_blobs];

        // A 3D matrix holding pointers to the products of the above maps to make an 'identity' map i.e.:
        //kinetic_double_map[blob_index][conf_index][via_conf_index] = kinetic_map[blob_index][via_conf_index][conf_index] * kinetic_map[blob_index][conf_index][via_conf_index]
        // Used to compare energies between conformers before actually switching
        kinetic_return_map = new SparseMatrixFixedPattern ***[params.num_blobs];

        // A 2D matrix holding the information about the 'num_states' kinetic states for each blob
        kinetic_state = new KineticState *[params.num_blobs];

        // Assign all memory
        for (i = 0; i < params.num_blobs; ++i)
        {
            kinetic_map[i] = new SparseMatrixFixedPattern *[params.num_conformations[i]];
            kinetic_return_map[i] = new SparseMatrixFixedPattern **[params.num_conformations[i]];

            if (params.num_conformations[i] == 1)
            {
                kinetic_map[i] = NULL;
                kinetic_return_map[i] = NULL;
            }
            else
            {
                for (j = 0; j < params.num_conformations[i]; ++j)
                {
                    kinetic_map[i][j] = new SparseMatrixFixedPattern[params.num_conformations[i]];
                    kinetic_return_map[i][j] = new SparseMatrixFixedPattern *[params.num_conformations[i]];
                }
            }
        }
    }

    // Load the vdw forcefield params matrix
    if (params.calc_ssint == 1 || params.calc_steric == 1)
    {
        if (ssint_matrix.init(params.ssint_in_fname, params.ssint_type, params.calc_ssint, params.ssint_cutoff) == FFEA_ERROR)
            FFEA_ERROR_MESSG("Error when reading from vdw forcefield params file.\n");
    }

    if (params.calc_vdw_rod == 1)
    {
        if (rod_lj_matrix.init(params.rod_lj_in_fname, "ljsteric", params.calc_vdw_rod, params.ssint_cutoff) == FFEA_ERROR)
            FFEA_ERROR_MESSG("Error when reading from rod Lennard-Jones params file (.rodlj).\n");

    }

    // detect how many threads we have for openmp
#ifdef USE_OPENMP
    num_threads = omp_get_max_threads();
    printf("\n\tNumber of threads detected: %d\n\n", num_threads);
#else
    num_threads = 1;
#endif

    // RNG: initialise the system, and allocate the Seeds:
    //    This has to be done before calling read_and_build_system.
    if (params.restart == 0)
    {
        // RNG.1 - allocate the Seeds:
        num_seeds = num_threads;
        if (params.calc_kinetics == 1)
            num_seeds += 1;
        Seeds = new unsigned long *[num_seeds];
        for (int i = 0; i < num_seeds; i++)
        {
            Seeds[i] = new unsigned long[6];
        }
        // RNG.2 - Initialise the package:
        // We need one rng for each thread, to avoid concurrency problems,
        // so generate an array of instances of RngStreams:
        // We first initialize the six seeds that RngStream needs:
        unsigned long sixseed[6];
        unsigned long twoMax[2] = {4294967087, 4294944443}; // these are max values for RngStream
        srand(params.rng_seed);
        for (int i = 0; i < 6; i++)
        {
            sixseed[i] = (rand() + rand()) % twoMax[i / 3];
        }
        // now initialise the package:
        RngStream::SetPackageSeed(sixseed);
        if (userInfo::verblevel > 1)
        {
            cout << "RngStream initialised using: ";
            for (int ni = 0; ni < 6; ni++)
            {
                cout << sixseed[ni] << " ";
            }
            cout << endl;
        }
        // RNG.3 - initialise the rngs related to the fluctuating stress, and noise.
        rng = new RngStream[num_threads];
        if (userInfo::verblevel > 2)
        {
            for (int ni = 0; ni < num_threads; ni++)
            {
                cout << "RNG[" << ni << "] initial state:" << endl;
                rng[ni].WriteState();
            }
        }
        // RNG.4 - and optionally initialise an extra Stream for kinetics.
        if (params.calc_kinetics == 1)
        {
            kinetic_rng = new RngStream;
            cout << "RngStream for kinetics created with initial state:" << endl;
            kinetic_rng->WriteState();
        }
    }
    else if (params.restart == 1)
    {
        // RNG - We'll now recover the state of the RNGs.
        printf("Getting state information from %s\n", params.icheckpoint_fname.c_str());
        // RNG.1 - READ Seeds FROM i.fcp into Seeds:
        // RNG.1.1 - open checkpoint file and check:
        ifstream checkpoint_in;
        checkpoint_in.open(params.icheckpoint_fname, ifstream::in);
        if (checkpoint_in.fail())
        {
            FFEA_FILE_ERROR_MESSG(params.icheckpoint_fname.c_str())
        }
        // RNG.1.2 - readlines, and num_seeds:
        vector<string> checkpoint_v;
        string line;
        while (getline(checkpoint_in, line))
        {
            checkpoint_v.push_back(line);
        }
        vector<string> header;
        int num_seeds_read;
        boost::split(header, checkpoint_v[0], boost::is_any_of(" "));
        try
        {
            num_seeds_read = stoi(header.back());
        }
        catch (invalid_argument &ia)
        {
            FFEA_ERROR_MESSG("Error reading the number of stress seeds: %s\n", ia.what());
        }
        int num_stress_seeds_read = num_seeds_read;
        if (params.calc_kinetics)
            num_seeds_read += 1;
        // RNG.1.3 - allocate Seeds:
        int num_active_rng = num_threads;
        if (params.calc_kinetics)
            num_active_rng += 1;
        int nlines = checkpoint_v.size();

        num_seeds = max(num_active_rng, num_seeds_read);

        // Allocate Seeds:
        Seeds = new unsigned long *[num_seeds];
        for (int i = 0; i < num_seeds; ++i)
        {
            Seeds[i] = new unsigned long[6];
            fill_n(Seeds[i], 6, 0); // fill the array with zeroes.
        }
        // RNG.1.4 - get Seeds :
        int cnt_seeds = 0;
        for (int i = 1; i < num_stress_seeds_read + 1; ++i)
        {
            vector<string> vline;
            boost::split(vline, checkpoint_v[i], boost::is_any_of(" "));
            // there must be 6 integers per line:
            if (vline.size() != 6)
            {
                FFEA_ERROR_MESSG("ERROR reading seeds\n")
            }
            for (int j = 0; j < 6; j++)
            {
                try
                {
                    Seeds[cnt_seeds][j] = stoul(vline[j]);
                }
                catch (invalid_argument &ia)
                {
                    FFEA_ERROR_MESSG("Error reading seeds as integers: %s\n", ia.what());
                }
            }
            cnt_seeds += 1;
        }

        // If kinetics active, one more seed to get (on line num_threads + 2)
        //   Trying to get this seed is essential: if we fail, it means that
        //   the previous run did not use kinetics, which is assumed later
        //   when initialising the RngStream for the kinetics.
        if (params.calc_kinetics)
        {
            vector<string> vline;
            boost::split(vline, checkpoint_v[num_threads + 2], boost::is_any_of(" "));
            // there must be 6 integers per line:
            if (vline.size() != 6)
            {
                FFEA_ERROR_MESSG("ERROR reading seeds\n")
            }
            for (int j = 0; j < 6; j++)
            {
                try
                {
                    Seeds[cnt_seeds][j] = stoul(vline[j]);
                }
                catch (invalid_argument &ia)
                {
                    FFEA_ERROR_MESSG("Error reading seeds as integers: %s\n", ia.what());
                }
            }
        }

        // RNG.2 - AND initialise rng:
        // RNG.2.1 - first the package and, and rng[0],
        cout << "To be initialised by: " << Seeds[0][0] << " " << Seeds[0][1] << " " << Seeds[0][2] << " " << Seeds[0][3] << " " << Seeds[0][4] << " " << Seeds[0][5] << endl;
        RngStream::SetPackageSeed(Seeds[0]);
        rng = new RngStream[num_threads];
        // RNG.2.2 - and now the rest of them: rng[1:]:
        for (int i = 1; i < (min(num_threads, num_stress_seeds_read)); i++)
        {
            rng[i].SetSeed(Seeds[i]);
        }
        if (userInfo::verblevel > 2)
        {
            for (int ni = 0; ni < num_threads; ni++)
            {
                cout << "RNG[" << ni << "] initial state:" << endl;
                rng[ni].WriteState();
            }
        }
        // if num_threads == num_threads in all the previous runs, it was that easy.
        // if num_threads < num_threads in a previous run, we have enough seeds, and we'll keep the
        //    extra seeds to be saved in ocheckpoint_fname.
        // if num_threads > num_threads in a previous run, we initialised the extra RNG properly
        //		when we set SetPackageSeed(Seeds[0]).
        // RNG.2.3 - if kinetics, initialise the last one!
        if (params.calc_kinetics == 1)
        {
            kinetic_rng = new RngStream;
            kinetic_rng->SetSeed(Seeds[num_seeds_read - 1]);
        }
    }

    // Build system of blobs, conformations, kinetics etc
    if (read_and_build_system(script_vector) != 0)
    {
        FFEA_error_text();
        cout << "System could not be built in World::read_and_build_system()" << endl;
        return FFEA_ERROR;
    }

    // If requested, initialise the PreComp_solver.
    //   Because beads need to be related to elements, it is much easier if
    //   it is done before moving the blobs to the latest trajectory step in
    //   case of "restart".
    if (params.calc_preComp == 1)
    {
        if (pc_solver.init(&pc_params, &params, blob_array) == FFEA_ERROR)
        {
            cout << "Failed to initialise PreComp_solver" << endl;
            return FFEA_ERROR;
        }
    }

    // Set up a Box:
    vector3 world_centroid, shift;
    get_system_centroid(&world_centroid);
    if (params.es_N_x < 1 || params.es_N_y < 1 || params.es_N_z < 1)
    {
        vector3 dimension_vector;
        get_system_dimensions(&dimension_vector);

        // Calculate decent box size
        params.es_N_x = 2 * (int)ceil(dimension_vector.x / params.ssint_cutoff);
        params.es_N_y = 2 * (int)ceil(dimension_vector.y / params.ssint_cutoff);
        params.es_N_z = 2 * (int)ceil(dimension_vector.z / params.ssint_cutoff);
    }

    // Move to box centre (if it is a new simulation! Otherwise trajectory will already have taken care of the move)
    box_dim.x = params.ssint_cutoff * params.es_N_x;
    box_dim.y = params.ssint_cutoff * params.es_N_y;
    box_dim.z = params.ssint_cutoff * params.es_N_z;

    shift.x = box_dim.x / 2.0 - world_centroid.x;
    shift.y = box_dim.y / 2.0 - world_centroid.y;
    shift.z = box_dim.z / 2.0 - world_centroid.z;

    std::cout << "Box dimensions: (" << box_dim.x << ", " << box_dim.y << ", " << box_dim.z << ")\n";
    std::cout << "World centroid: (" << world_centroid.x << ", " << world_centroid.y << ", " << world_centroid.z << ")\n";
    std::cout << "Shift: (" << shift.x << ", " << shift.y << ", " << shift.z << ")\n";

    if (params.move_into_box == 1)
    { // && params.restart == 0)
        for (i = 0; i < params.num_blobs; i++)
        {
            //active_blob_array[i]->get_centroid(&world_centroid);
            active_blob_array[i]->move(shift.x, shift.y, shift.z);
            active_blob_array[i]->calc_all_centroids();
        }
        float shift_rod[3] = {(float)shift.x, (float)shift.y, (float)shift.z}; // this class is some historical junk
        for (i = 0; i < params.num_rods; i++)
        {
            if (rod_array[i]->restarting == false)
            {
                rod_array[i]->translate_rod(rod_array[i]->current_r, shift_rod);
                rod_array[i]->translate_rod(rod_array[i]->equil_r, shift_rod);
                rod_array[i]->write_frame_to_file(); // rod traj contains initial state but positioned in box
            }
        }
        // todo: for each attachment, move the attachment_node_pos by shift
        for (i = 0; i < params.num_interfaces; i++)
        {
            //for(int j=0; j<3; j++){rod_blob_interface_array[i]->attachment_node_pos[j] += shift_rod[j];}
            rod_blob_interface_array[i]->update_internal_state(true, true);
            //            if(rod::dbg_print){rod::print_array("shifted interface position", rod_blob_interface_array[i]->attachment_node_pos_equil, 3);}
        }

        if (rod::dbg_print)
        {
            rod::print_array("shift", shift_rod, 3);
        }
    }

    // Check that rods aren't too big
    if (params.pbc_rod == 1)
    {
        for (int i = 0; i < params.num_rods; i++)
            rod_box_length_check(rod_array[i], {(float)box_dim.x, (float)box_dim.y, (float)box_dim.z});
    }

    // Now everything has been moved into boxes etc, save all initial positions
    for (i = 0; i < params.num_blobs; ++i)
    {
        for (j = 0; j < params.num_conformations[i]; ++j)
        {
            blob_array[i][j].set_pos_0();
        }
    }
    for (int i = 0; i < params.num_interfaces; i++)
    {
        rod_blob_interface_array[i]->update_J_0();
    }

    // If not restarting a previous simulation, create new trajectory and measurement files. But only if full simulation is happening!
    if (mode == 0)
    {
        // In any case, open the output checkpoint file for writing
        if ((checkpoint_out = fopen(params.ocheckpoint_fname.c_str(), "wb")) == NULL)
        {
            FFEA_FILE_ERROR_MESSG(params.ocheckpoint_fname.c_str());
        }
#ifdef FFEA_PARALLEL_FUTURE
        // And launch a first trajectory thread, that will be catched up at print_traj time
        thread_writingTraj = std::async(std::launch::async, &World::do_nothing, this);
#endif

        if (params.restart == 0)
        {

            // Open the trajectory output file for writing
            if ((trajectory_out = fopen(params.trajectory_out_fname.c_str(), "wb")) == NULL)
            {
                FFEA_FILE_ERROR_MESSG(params.trajectory_out_fname.c_str())
            }

            // HEADER FOR TRAJECTORY
            // Print initial info stuff
            fprintf(trajectory_out, "FFEA_trajectory_file\n\nInitialisation:\nNumber of Blobs %d\nNumber of Conformations", params.num_blobs);
            for (i = 0; i < params.num_blobs; ++i)
            {
                fprintf(trajectory_out, " %d", params.num_conformations[i]);
            }
            fprintf(trajectory_out, "\n");

            for (i = 0; i < params.num_blobs; ++i)
            {
                fprintf(trajectory_out, "Blob %d:\t", i);
                for (j = 0; j < params.num_conformations[i]; ++j)
                {
                    fprintf(trajectory_out, "Conformation %d Nodes %d\t", j, blob_array[i][j].get_num_nodes());
                }
                fprintf(trajectory_out, "\n");
            }
            fprintf(trajectory_out, "\n");

            // First line in trajectory data should be an asterisk (used to delimit different steps for easy seek-search in restart code)
            fprintf(trajectory_out, "*\n");

            // Open the measurement output file for writing
            if ((measurement_out = fopen(params.measurement_out_fname.c_str(), "wb")) == NULL)
            {
                FFEA_FILE_ERROR_MESSG(params.measurement_out_fname.c_str())
            }

            // HEADER FOR MEASUREMENTS
            // Write header to output file
            write_output_header(measurement_out, FFEA_script_filename);

            // Write params to this output file
            params.write_to_file(measurement_out, pc_params);

            // Get ready to write the measurements (this is the order things must be written later. There will be no floating zeroes!)
            fprintf(measurement_out, "Measurements:\n");
            fprintf(measurement_out, "%-14s", "Time");

            // Do we need kinetic energy?
            if (mass_in_system)
            {
                fprintf(measurement_out, "%-14s", "KineticEnergy");
            }
            fprintf(measurement_out, "%-14s", "StrainEnergy");
            fprintf(measurement_out, "%-14s%-14s%-14s%-14s", "Centroid.x", "Centroid.y", "Centroid.z", "RMSD");

            // Are these field enabled?
            if (params.calc_springs != 0)
            {
                fprintf(measurement_out, "%-14s", "SpringEnergy");
            }
            if (params.calc_ssint == 1 || params.calc_steric == 1)
            {
                fprintf(measurement_out, "%-15s", "SurfSurfEnergy");
            }
            if (params.calc_preComp != 0)
            {
                fprintf(measurement_out, "%-14s", "PreCompEnergy");
            }
            fprintf(measurement_out, "\n");
            fflush(measurement_out);

            // HEADER FOR DETAILED MEASUREMENTS (if necessary)
            if (writeDetailed)
            {
                detailed_meas_out = fopen(params.detailed_meas_out_fname.c_str(), "wb");
                fprintf(detailed_meas_out, "FFEA Detailed Measurement File\n\nMeasurements:\n");
                fprintf(detailed_meas_out, "%-14s", "Time");
                for (i = 0; i < params.num_blobs; ++i)
                {
                    fprintf(detailed_meas_out, "| B%d ", i);
                    if (active_blob_array[i]->there_is_mass())
                    {
                        fprintf(detailed_meas_out, "%-14s", "KineticEnergy");
                    }
                    fprintf(detailed_meas_out, "%-14s", "StrainEnergy");
                    fprintf(detailed_meas_out, "%-14s%-14s%-14s%-14s", "Centroid.x", "Centroid.y", "Centroid.z", "RMSD");
                }

                if (params.calc_ssint == 1 || params.calc_steric == 1 || params.calc_preComp == 1 || params.calc_springs == 1)
                {
                    for (i = 0; i < params.num_blobs; ++i)
                    {
                        for (j = i; j < params.num_blobs; ++j)
                        {
                            fprintf(detailed_meas_out, "| B%dB%d ", i, j);
                            if (active_blob_array[i]->there_is_ssint() && active_blob_array[j]->there_is_ssint())
                            {
                                fprintf(detailed_meas_out, "%-15s", "SurfSurfEnergy");
                            }
                            if (active_blob_array[i]->there_are_springs() && active_blob_array[j]->there_are_springs())
                            {
                                fprintf(detailed_meas_out, "%-14s", "SpringEnergy");
                            }

                            if (active_blob_array[i]->there_are_beads() && active_blob_array[j]->there_are_beads())
                            {
                                fprintf(detailed_meas_out, "%-14s", "PreCompEnergy");
                            }
                        }
                    }
                }
                fprintf(detailed_meas_out, "\n");
                fflush(detailed_meas_out);
            }

            // Open the beads_out file:
            if (params.trajbeads_fname_set == 1)
            {
                if ((trajbeads_out = fopen(params.trajectory_beads_fname.c_str(), "wb")) == NULL)
                {
                    FFEA_FILE_ERROR_MESSG(params.trajectory_beads_fname.c_str())
                }
            }

            // Open the kinetics output file for writing (if neccessary) and write initial stuff
            if (params.kinetics_out_fname_set == 1)
            {
                if ((kinetics_out = fopen(params.kinetics_out_fname.c_str(), "wb")) == NULL)
                {
                    FFEA_FILE_ERROR_MESSG(params.kinetics_out_fname.c_str())
                }
                fprintf(kinetics_out, "FFEA_kinetic_trajectory_file\n\nNumber of Blobs %d\n\n", params.num_blobs);
                for (i = 0; i < params.num_blobs; ++i)
                {
                    fprintf(kinetics_out, "              Blob %d         ", i);
                }
                fprintf(kinetics_out, "\n# step ");
                for (i = 0; i < params.num_blobs; ++i)
                {
                    fprintf(kinetics_out, "|state | conformation |");
                }
                fprintf(kinetics_out, "\n");
                fflush(kinetics_out);
            }
        }
        else
        {

            // Otherwise, seek backwards from the end of the trajectory file looking for '*' character (delimitter for snapshots)

            /*
             * Trajectory first
             */
            bool singleframe = false;
            char c;

            printf("Restarting from trajectory file %s\n", params.trajectory_out_fname.c_str());
            if ((trajectory_out = fopen(params.trajectory_out_fname.c_str(), "rb")) == NULL)
            {
                FFEA_FILE_ERROR_MESSG(params.trajectory_out_fname.c_str())
            }

            printf("Reverse searching for 4 asterisks ");
            if (frames_to_delete != 0)
            {
                printf(", plus an extra %d, ", (frames_to_delete)*2);
            }
            printf("(denoting %d completely written snapshots)...\n", frames_to_delete + 1);
            if (fseek(trajectory_out, 0, SEEK_END) != 0)
            {
                FFEA_ERROR_MESSG("Could not seek to end of file\n")
            }

            // Variable to store position of last asterisk in trajectory file (initialise it at end of file)
            long last_asterisk_pos;
            last_asterisk_pos = ftell(trajectory_out);

            int num_asterisks = 0;
            int num_asterisks_to_find = 3 + (frames_to_delete)*2 + 1; // 3 to get to top of last frame, then two for every subsequent frame. Final 1 to find ending conformations of last step
            while (num_asterisks != num_asterisks_to_find)
            {
                if (fseek(trajectory_out, -2, SEEK_CUR) != 0)
                {
                    //perror(NULL);
                    //FFEA_ERROR_MESSG("It is likely we have reached the begininng of the file whilst trying to delete frames. You can't delete %d frames.\n", frames_to_delete)
                    printf("Found beginning of file. Searching forwards for next asterisk...");
                    singleframe = true;

                    // This loop will allow the script to find the 'final' asterisk
                    while (true)
                    {
                        if ((c = fgetc(trajectory_out)) == '*')
                        {
                            fseek(trajectory_out, -1, SEEK_CUR);
                            break;
                        }
                    }
                }
                c = fgetc(trajectory_out);
                if (c == '*')
                {
                    num_asterisks++;
                    printf("Found %d\n", num_asterisks);

                    // get the position in the file of this last asterisk
                    //if (num_asterisks == num_asterisks_to_find - 2) {
                    //  last_asterisk_pos = ftello(trajectory_out);
                    //}
                }
            }

            // char sline[255];
            if ((c = fgetc(trajectory_out)) != '\n')
            {
                ungetc(c, trajectory_out);
            }
            else
            {
                last_asterisk_pos = ftell(trajectory_out);
            }

            // Get the conformations for the last snapshot (or set them as 0 if we have only 1 frame)
            int current_conf, crap;
            if (!singleframe)
            {
                crap = fscanf(trajectory_out, "Conformation Changes:\n");
                for (i = 0; i < params.num_blobs; ++i)
                {
                    crap = fscanf(trajectory_out, "Blob %*d: Conformation %*d -> Conformation %d\n", &current_conf);
                    active_blob_array[i] = &blob_array[i][current_conf];
                }
                crap = fscanf(trajectory_out, "*\n");
                last_asterisk_pos = ftell(trajectory_out);
            }
            else
            {
                for (i = 0; i < params.num_blobs; ++i)
                {
                    active_blob_array[i] = &blob_array[i][0];
                }
            }

            // Load this frame
            printf("Loading Blob position and velocity data from last completely written snapshot \n");
            int blob_id, conformation_id;
            long long rstep;
            for (int b = 0; b < params.num_blobs; b++)
            {
                if (fscanf(trajectory_out, "Blob %d, Conformation %d, step %lld\n", &blob_id, &conformation_id, &rstep) != 3)
                {
                    FFEA_ERROR_MESSG("Error reading header info for Blob %d\n", b)
                }
                if (blob_id != b)
                {
                    FFEA_ERROR_MESSG("Mismatch in trajectory file - found blob id = %d, was expecting blob id = %d\n", blob_id, b)
                }
                printf("Loading node position, velocity and potential from restart trajectory file, for blob %d, step %lld\n", blob_id, rstep);
                if (active_blob_array[b]->read_nodes_from_file(trajectory_out) == FFEA_ERROR)
                {
                    FFEA_ERROR_MESSG("Error restarting blob %d\n", blob_id)
                }
            }

            // Final conformation bit
            crap = fscanf(trajectory_out, "*\nConformation Changes:\n");
            for (i = 0; i < params.num_blobs; ++i)
            {
                crap = fscanf(trajectory_out, "Blob %*d: Conformation %*d -> Conformation %*d\n");
            }
            crap = fscanf(trajectory_out, "*\n");

            // Set truncation location
            //   last_asterisk_pos = ftell(trajectory_out);
            step_initial = rstep;
            printf("...done. Simulation will commence from step %lld\n", step_initial);
            fclose(trajectory_out);

            // Truncate the trajectory file up to the point of the stored asterisk position, thereby deleting the last frame and erasing any half-written time steps that may occur after it
            printf("Truncating the trajectory file to the last asterisk...\n");
            error_code ec1;
            filesystem::resize_file(params.trajectory_out_fname, last_asterisk_pos, ec1);
            if (ec1)
            {
                FFEA_ERROR_MESSG("Error when trying to truncate trajectory file %s\n", params.trajectory_out_fname.c_str())
            }

            /*
             * Measurement files
             */

            // Global

            if ((measurement_out = fopen(params.measurement_out_fname.c_str(), "rb")) == NULL)
            {
                FFEA_FILE_ERROR_MESSG(params.measurement_out_fname.c_str())
            }

            if (fseek(measurement_out, 0, SEEK_END) != 0)
            {
                FFEA_ERROR_MESSG("Could not seek to end of file\n")
            }

            last_asterisk_pos = ftell(measurement_out);

            // Looking for newlines this time, as each measurement frame is a single line
            int num_newlines = 0;
            int num_newlines_to_find = frames_to_delete + 1; // 1 for every frame only, and 1 extra as we redo the last step. No need to read in the last meas line
            while (num_newlines != num_newlines_to_find)
            {
                if (fseek(measurement_out, -2, SEEK_CUR) != 0)
                {
                    FFEA_ERROR_MESSG("Error when trying to find last frame from file %s\n", params.measurement_out_fname.c_str())
                }
                c = fgetc(measurement_out);
                if (c == '\n')
                {
                    num_newlines++;
                    printf("Found %d\n", num_newlines);
                }
            }

            last_asterisk_pos = ftell(measurement_out);

            // Truncate the measurement file up to the point of the last newline
            printf("Truncating the measurement file to the appropriate line...\n");
            error_code ec2;
            filesystem::resize_file(params.measurement_out_fname, last_asterisk_pos, ec2);
            if (ec2)
            {
                FFEA_ERROR_MESSG("Error when trying to truncate measurment file %s\n", params.measurement_out_fname.c_str())
            }

            // Append a newline to the end of this truncated measurement file (to replace the one that may or may not have been there)
            //fprintf(measurement_out, "#==RESTART==\n");

            last_asterisk_pos = ftell(measurement_out);

            // Detailed
            if (writeDetailed)
            {

                if ((detailed_meas_out = fopen(params.detailed_meas_out_fname.c_str(), "rb")) == NULL)
                {
                    FFEA_FILE_ERROR_MESSG(params.detailed_meas_out_fname.c_str())
                }

                if (fseek(detailed_meas_out, 0, SEEK_END) != 0)
                {
                    FFEA_ERROR_MESSG("Could not seek to end of file\n")
                }

                last_asterisk_pos = ftell(detailed_meas_out);

                // Looking for newlines this time, as each measurement frame is a single line
                num_newlines = 0;
                num_newlines_to_find = frames_to_delete + 1; // 1 for every frame, plus the first one, assuming all were written correctly
                while (num_newlines != num_newlines_to_find)
                {
                    if (fseek(detailed_meas_out, -2, SEEK_CUR) != 0)
                    {
                        FFEA_ERROR_MESSG("Error when trying to find last frame from file %s\n", params.detailed_meas_out_fname.c_str())
                    }
                    c = fgetc(detailed_meas_out);
                    if (c == '\n')
                    {
                        num_newlines++;
                        printf("Found %d\n", num_newlines);
                    }
                }

                last_asterisk_pos = ftell(detailed_meas_out);

                // Truncate the measurement file up to the point of the last newline
                printf("Truncating the detailed measurement file to the appropriate line...\n");
                error_code ec3;
                filesystem::resize_file(params.detailed_meas_out_fname, last_asterisk_pos, ec3);
                if (ec3)
                {
                    FFEA_ERROR_MESSG("Error when trying to truncate measurment file %s\n", params.detailed_meas_out_fname.c_str())
                }
            }

            // Open trajectory and measurment files for appending
            printf("Opening trajectory and measurement files for appending.\n");
            if ((trajectory_out = fopen(params.trajectory_out_fname.c_str(), "ab")) == NULL)
            {
                FFEA_FILE_ERROR_MESSG(params.trajectory_out_fname.c_str())
            }
            if ((measurement_out = fopen(params.measurement_out_fname.c_str(), "ab")) == NULL)
            {
                FFEA_FILE_ERROR_MESSG(params.measurement_out_fname.c_str())
            }

            // And the detailed meas file, maybe
            if (writeDetailed)
            {
                if ((detailed_meas_out = fopen(params.detailed_meas_out_fname.c_str(), "ab")) == NULL)
                {
                    FFEA_FILE_ERROR_MESSG(params.detailed_meas_out_fname.c_str())
                }
            }

            // And kinetic file
            if (params.kinetics_out_fname_set == 1)
            {
                if ((kinetics_out = fopen(params.kinetics_out_fname.c_str(), "ab")) == NULL)
                {
                    FFEA_FILE_ERROR_MESSG(params.kinetics_out_fname.c_str())
                }
                //fprintf(kinetics_out, "#==RESTART==\n");
            }

            /*
            *
            *
            * Fix restart for measurements and kinetics in future. Use rstep to find appropriate line
            *
            */
            // traj beads do not work yet:
            if (params.trajbeads_fname_set == 1)
                FFEA_ERROR_MESSG("FFEA cannot still restart and keep writing on the beads file. Just remove it from your input file.");
        }
    }

    // Check if there are static blobs:
    bool there_are_static_blobs = false;
    for (i = 0; i < params.num_blobs; i++)
    {
        for (j = 0; j < params.num_conformations[i]; ++j)
        {
            if (blob_array[i][j].get_motion_state() != FFEA_BLOB_IS_DYNAMIC)
            {
                there_are_static_blobs = true;
                break;
            }
        }
        if (there_are_static_blobs == true)
            break;
    }

    // 'calc_steric' and 'ssint_type' are now independent. Can the solvers be refactored to reduce on code?
    if (params.calc_ssint == 1 || params.calc_steric == 1)
    {
        if (params.ssint_type == "lennard-jones")
            vdw_solver = new (std::nothrow) VdW_solver();
        else if (params.ssint_type == "steric")
            vdw_solver = new (std::nothrow) Steric_solver();
        else if (params.ssint_type == "ljsteric")
            vdw_solver = new (std::nothrow) LJSteric_solver();
        else if (params.ssint_type == "gensoft")
            vdw_solver = new (std::nothrow) GenSoftSSINT_solver();

        if (vdw_solver == NULL)
            FFEA_ERROR_MESSG("World::init failed to initialise the ssint_solver.\n");

        vdw_solver->init(&lookup, &box_dim, &ssint_matrix, params.steric_factor, params.num_blobs, params.inc_self_ssint, params.ssint_type, params.steric_dr, params.calc_kinetics, there_are_static_blobs);
    }

    // Calculate the total number of vdw interacting faces in the entire system
    total_num_surface_faces = 0;
    for (i = 0; i < params.num_blobs; i++)
    {
        for (j = 0; j < params.num_conformations[i]; ++j)
        {
            total_num_surface_faces += blob_array[i][j].get_num_faces();
        }
    }
    printf("Total number of surface faces in system: %d\n", total_num_surface_faces);

    // Initialise the face-face neighbour list for vdw or es:
    if (params.calc_ssint == 1 || params.calc_steric == 1 || params.calc_es == 1)
    {

        // Allocate memory for an NxNxN grid, with a 'pool' for the required number of surface faces
        printf("Allocating memory for nearest neighbour lookup grid...\n");
#ifdef FFEA_PARALLEL_FUTURE
        int lookup_error = lookup.alloc_dual(params.es_N_x, params.es_N_y, params.es_N_z, total_num_surface_faces);
#else
        int lookup_error = lookup.alloc(params.es_N_x, params.es_N_y, params.es_N_z, total_num_surface_faces);
#endif
        if (lookup_error == FFEA_ERROR)
        {
            FFEA_error_text();
            printf("When allocating memory for nearest neighbour lookup grid\n");
            return FFEA_ERROR;
        }
        printf("...done\n");

        printf("Box has volume %e cubic angstroms\n", (box_dim.x * box_dim.y * box_dim.z) * mesoDimensions::volume * 1e30);

        // Add all the faces from each Blob to the lookup pool
        printf("Adding all faces to nearest neighbour grid lookup pool\n");
        for (i = 0; i < params.num_blobs; i++)
        {
            for (j = 0; j < params.num_conformations[i]; ++j)
            {
                int num_faces_added = 0;
                for (k = 0; k < blob_array[i][j].get_num_faces(); k++)
                {
                    Face *b_face = blob_array[i][j].get_face(k);
                    if (b_face != NULL)
                    {
#ifdef FFEA_PARALLEL_FUTURE
                        lookup_error = lookup.add_to_pool_dual(b_face);
#else
                        lookup_error = lookup.add_to_pool(b_face);
#endif
                        if (lookup_error == FFEA_ERROR)
                        {
                            FFEA_error_text();
                            printf("When attempting to add a face to the lookup pool\n");
                            return FFEA_ERROR;
                        }
                        num_faces_added++;
                    }
                }
                if (userInfo::verblevel > 1)
                    printf("%d 'ssint active' faces, from blob %d, conformation %d, added to lookup grid.\n", num_faces_added, i, j);
            }
        }

#ifdef FFEA_PARALLEL_FUTURE
        // And build the lookup table for the first time:
        thread_updatingVdWLL = std::async(std::launch::async, &World::prebuild_nearest_neighbour_lookup_wrapper, this, params.ssint_cutoff);
#endif
    }

    // pc_solver has already allocated memory for  its own neighbour list.
    // Still it had to wait until everything was put into box,
    //   to place the beads onto the voxels.
    // We now calculate the bead positions
    //    to be able to build_pc_nearest_neighbour_lookup()
    if (params.calc_preComp == 1)
    {
        pc_solver.compute_bead_positions();
#ifdef FFEA_PARALLEL_FUTURE
        thread_updatingPCLL = std::async(std::launch::async, &PreComp_solver::prebuild_pc_nearest_neighbour_lookup, &pc_solver);
#endif
    }

    // Initialise the BEM PBE solver
    if (params.calc_es == 1)
    {
        printf("Initialising Boundary Element Poisson Boltzmann Solver...\n");
        PB_solver.init(&lookup);
        PB_solver.set_kappa(params.kappa);

        // Initialise the nonsymmetric matrix solver
        nonsymmetric_solver.init(total_num_surface_faces, params.epsilon2, params.max_iterations_cg);

        // Allocate memory for surface potential vector
        phi_Gamma = new scalar[total_num_surface_faces];
        for (i = 0; i < total_num_surface_faces; i++)
            phi_Gamma[i] = 0;

        work_vec = new scalar[total_num_surface_faces];
        for (i = 0; i < total_num_surface_faces; i++)
            work_vec[i] = 0;

        J_Gamma = new scalar[total_num_surface_faces];
        for (i = 0; i < total_num_surface_faces; i++)
            J_Gamma[i] = 0;
    }

#ifdef FFEA_PARALLEL_WITHIN_BLOB
    printf("Now initialised with 'within-blob parallelisation' (FFEA_PARALLEL_WITHIN_BLOB) on %d threads.\n", num_threads);
#endif

#ifdef FFEA_PARALLEL_PER_BLOB
    printf("Now initialised with 'per-blob parallelisation' (FFEA_PARALLEL_PER_BLOB) on %d threads.\n", num_threads);
#endif

#ifdef USE_MPI
    et = MPI::Wtime() - st;
    cout << "benchmarking--------Initialising time of ffea :" << et << "seconds" << endl;
#endif

    return FFEA_OK;
}

/**
 * @brief Finds the largest allowed timesteps
 * @details By linearising the equation of motion, this function performs matrix
 *   algebra using the Eigen libraries to find the largest allowed timestep for
 *   ffea numerical integration.
 * */
int World::get_smallest_time_constants()
{

    // This is currently only for active blob, as inactive blobs are all at infinity due to linked list problems

    // Global variables
    int i, j, dt_min_bin, dt_max_bin, num_nodes, num_rows;
    scalar dt_min_world = INFINITY, dt_max_world = -1 * INFINITY;
    string dt_min_world_type = "viscous", dt_max_world_type = "viscous";
    Eigen::EigenSolver<Eigen_MatrixX> es_v, es_m;
    vector<scalar> tauv, taum;
    vector<scalar>::iterator it;

    cout << "Calculating time constants..." << endl
         << endl;
    for (i = 0; i < params.num_blobs; ++i)
    {
        cout << "\tBlob " << i << ":" << endl
             << endl;

        // Ignore if we have a static blob
        if (active_blob_array[i]->get_motion_state() == FFEA_BLOB_IS_STATIC)
        {
            cout << "\t\tBlob " << i << " is STATIC. No associated timesteps." << endl;
            continue;
        }

        // What will be the fastest dynamics? Inertial or viscous?

        // Viscous only

        // Build matrices
        num_nodes = active_blob_array[i]->get_num_linear_nodes();

        // Direction matters here
        num_rows = 3 * num_nodes;

        Eigen::SparseMatrix<scalar> K(num_rows, num_rows);
        Eigen::SparseMatrix<scalar> A(num_rows, num_rows);
        Eigen_MatrixX K_inv(num_rows, num_rows);
        Eigen_MatrixX I(num_rows, num_rows);
        Eigen_MatrixX tau_inv(num_rows, num_rows);

        // Build viscosity matrix, K
        cout << "\r\t\tCalculating the Viscosity Matrix, K (task 1/5)..." << flush;
        if (active_blob_array[i]->build_linear_node_viscosity_matrix(&K) == FFEA_ERROR)
        {
            cout << endl
                 << "\t\t";
            FFEA_error_text();
            cout << endl
                 << endl
                 << "In function 'Blob::build_linear_node_viscosity_matrix' from blob " << i << endl;
            return FFEA_ERROR;
        }
        cout << "done!" << flush;
        //cout << K * (mesoDimensions::force / mesoDimensions::velocity) << endl;
        // Build elasticity matrix, A
        cout << "\r\t\tCalculating the Elasticity Matrix, A (task 2/5)..." << flush;
        if (active_blob_array[i]->build_linear_node_elasticity_matrix(&A) == FFEA_ERROR)
        {
            cout << endl
                 << "\t\t";
            FFEA_error_text();
            cout << "In function 'Blob::build_linear_node_elasticity_matrix'" << i << endl;
            return FFEA_ERROR;
        }
        cout << "done!" << flush;

        // Invert K (it's symmetric! Will not work if stokes_visc == 0)
        cout << "\r\t\tAttempting to invert K to form K_inv (task 3/5)..." << flush;
        Eigen::SimplicialCholesky<Eigen::SparseMatrix<scalar>> Cholesky(K); // performs a Cholesky factorization of K
        I.setIdentity();
        K_inv = Cholesky.solve(I);
        if (Cholesky.info() == Eigen::Success)
        {
            cout << "done!" << flush;
        }
        else if (Cholesky.info() == Eigen::NumericalIssue)
        {
            cout << endl
                 << "\t\t"
                 << "Viscosity Matrix could not be inverted via Cholesky factorisation due to numerical issues. You possibly don't have an external solvent set, or it is too low." << endl;
            return FFEA_ERROR;
        }
        else if (Cholesky.info() == Eigen::NoConvergence)
        {
            cout << "\nInversion iteration couldn't converge. K must be a crazy matrix. Possibly has zero eigenvalues?" << endl;
            return FFEA_ERROR;
        }

        // Apply to A
        cout << "\r\t\tCalculating inverse time constant matrix, tau_inv = K_inv * A (task 4/5)..." << flush;
        tau_inv = K_inv * A;
        cout << "done!" << flush;
        printf("%c[2K", 27);

        // Diagonalise
        //cout << "\r" << "\t\t                                                                     " << flush;
        cout << "\r\t\tDiagonalising tau_inv (task 5/5)..." << flush;
        es_v.compute(tau_inv);
        for (j = 0; j < num_rows; ++j)
        {
            tauv.push_back(1.0 / fabs(es_v.eigenvalues()[j].real()));
        }
        cout << "done!" << endl
             << flush;

        if (active_blob_array[i]->get_linear_solver() != FFEA_NOMASS_CG_SOLVER)
        {

            // Inertial 'always' fastest

            // Build matrices
            num_nodes = active_blob_array[i]->get_num_linear_nodes();

            // Direction still matters here due to viscosity (potentially)
            num_rows = 3 * num_nodes;
            Eigen::SparseMatrix<scalar> M(num_rows, num_rows);
            Eigen_MatrixX M_inv(num_rows, num_rows);
            Eigen::VectorXd Mev(num_rows);
            //Mev.setZero();
            //for(j = 0; j < num_nodes; ++j) {
            //	Mev(3*j) = 1.0;
            //}

            // Build mass matrix, M
            cout << "\r\t\tCalculating the Mass Matrix, M (task 1/4)..." << flush;
            if (active_blob_array[i]->build_linear_node_mass_matrix(&M) == FFEA_ERROR)
            {
                cout << endl
                     << "\t\t";
                FFEA_error_text();
                cout << endl
                     << endl
                     << "In function 'Blob::build_linear_node_viscosity_matrix' from blob " << i << endl;
                return FFEA_ERROR;
            }
            cout << "done!" << flush;
            //cout << endl << flush << mesoDimensions::mass * (Mev.transpose() * M * Mev) << endl;
            //cout << (4.0 / 3.0) * 3.1415926535 * 125e-27 * 1500 << endl;
            //exit(0);

            // Build eigenvector, apply and output the mass...

            // Invert M (it's symmetric!)
            cout << "\r\t\tAttempting to invert M to form M_inv (task 2/4)..." << flush;
            Eigen::SimplicialCholesky<Eigen::SparseMatrix<scalar>> Cholesky(M); // performs a Cholesky factorization of M
            M_inv = Cholesky.solve(I);
            if (Cholesky.info() == Eigen::Success)
            {
                cout << "done!" << flush;
            }
            else if (Cholesky.info() == Eigen::NumericalIssue)
            {
                cout << endl
                     << "\t\t"
                     << "Mass Matrix could not be inverted via Cholesky factorisation due to numerical issues. This...should not be the case. You have a very odd mass distribution. Try the CG_nomass solver" << endl;
                return FFEA_ERROR;
            }
            else if (Cholesky.info() == Eigen::NoConvergence)
            {
                cout << endl
                     << "\t\t"
                     << "Inversion iteration couldn't converge. M must be a crazy matrix. Possibly has zero eigenvalues? Try the CG_nomass solver." << endl;
                return FFEA_ERROR;
            }

            // Apply to K
            cout << "\r\t\tCalculating inverse time constant matrix, tau_inv = M_inv * K (task 3/4)..." << flush;
            tau_inv = M_inv * K;
            cout << "done!" << flush;
            printf("%c[2K", 27);
            // Diagonalise
            //cout << "\r" << "\t\t                                                                                                            " << flush;
            cout << "\r\t\tDiagonalising tau_inv (task 4/4)..." << flush;
            es_m.compute(tau_inv);
            for (j = 0; j < num_rows; ++j)
            {
                taum.push_back(1.0 / fabs(es_m.eigenvalues()[j].real()));
            }
            cout << "done!" << endl
                 << endl
                 << flush;
        }

        // But is it a numerical instability problem, or a small elements problem? Solve 1 step to find out (at a later date)

        // Get extreme timesteps from this eigendecomposition.)

        // Sort eigenvalues
        cout << "\r\t\tSorting eigenvalues (task 1/1)..." << flush;
        sort(tauv.begin(), tauv.end());
        sort(taum.begin(), taum.end());
        cout << "done!" << flush;

        // Ignore the 6 translational / rotational modes (they are likely the slowest 6 modes)
        scalar dt_max_blob = tauv.at(num_rows - 7);
        scalar dt_min_blob = tauv.at(0);
        string dt_min_blob_type = "viscous", dt_max_blob_type = "viscous";

        //for(int l = 0; l < 10; ++l) {
        //	cout << endl << "ev " << l << " visc - " << tauv.at(l) << " mass - " << taum.at(l) << endl;
        //}
        //for(int l = num_rows - 10; l < num_rows; ++l) {
        //		cout << endl << "ev " << l << " visc - " << tauv.at(l) << " mass - " << taum.at(l) << endl;
        //}
        //exit(0);

        // We don't need to ignore top 6 here as they have energy associated with them now i.e. not zero eigenvalues
        if (active_blob_array[i]->get_linear_solver() != FFEA_NOMASS_CG_SOLVER)
        {

            if (taum.at(num_rows - 1) > dt_max_blob)
            {
                dt_max_blob = taum.at(num_rows - 1);
                dt_max_blob_type = "inertial";
            }
            if (taum.at(0) < dt_min_blob)
            {
                dt_min_blob = taum.at(0);
                dt_min_blob_type = "inertial";
            }
        }

        //cout << "\r\t\tThe time-constant of the slowest mode in Blob " << blob_index << ", tau_max = " << (1.0 / smallest_val) * mesoDimensions::time << "s" << endl;
        //cout << "\t\tThe time-constant of the fastest mode in Blob " << blob_index << ", tau_min = " << (1.0 / largest_val) * mesoDimensions::time << "s" << endl << endl;
        cout << "\r\t\tFastest Mode: tau (" << dt_min_blob_type << ") = " << dt_min_blob * mesoDimensions::time << "s" << endl;
        cout << "\t\tSlowest Mode: tau (" << dt_max_blob_type << ") = " << dt_max_blob * mesoDimensions::time << "s" << endl
             << endl;

        // Global stuff
        if (dt_max_blob > dt_max_world)
        {
            dt_max_world = dt_max_blob;
            dt_max_world_type = dt_max_blob_type;
            dt_max_bin = i;
        }

        if (dt_min_blob < dt_min_world)
        {
            dt_min_world = dt_min_blob;
            dt_min_world_type = dt_min_blob_type;
            dt_min_bin = i;
        }
    }

    cout << endl
         << "Global Time Constant Details:" << endl
         << endl;
    cout << "\t\tFastest Mode: Blob " << dt_min_bin << ", tau (" << dt_min_world_type << ") = " << dt_min_world * mesoDimensions::time << "s" << endl;
    cout << "\t\tSlowest Mode: Blob " << dt_max_bin << ", tau (" << dt_max_world_type << ") = " << dt_max_world * mesoDimensions::time << "s" << endl
         << endl;
    if (dt_min_world_type == "inertial")
    {
        cout << "\t\tPlease make sure your simulation timestep is less than " << dt_min_world * mesoDimensions::time << "s, for a stable simulation." << endl;
        cout << "\t\tTake note than the energies will become inaccurate before this, so check your energy equilibrates correctly. If unsure, set dt << " << dt_min_world * mesoDimensions::time << "s" << endl
             << endl;
    }
    cout << "\t\tFor dynamical convergence, your simulation must run for longer than " << dt_max_world * mesoDimensions::time << "s." << endl
         << endl;

    cout << "\t\tFINAL NOTE - If, after taking into account the above time constants, your simulation still fails (due to element inversion) it is not due to numerical instability from the integration, ";
    cout << "but because a single timestep, with the average size of the noise, causes a step size larger than your smallest element. You must therefore coarsen your mesh further for the ";
    cout << "continuum approximation to be valid. Thanks :)" << endl
         << endl;
    return FFEA_OK;
}

/*int World::get_smallest_time_constants() {

	int blob_index;
	int num_nodes, num_rows;
	scalar dt_min = INFINITY;
	scalar dt_max = -1 * INFINITY;
	int dt_max_bin = -1, dt_min_bin = -1;
	cout << "\n\nFFEA mode - Calculate Maximum Allowed Timesteps" << endl << endl;

	// Do active blobs first
	for(blob_index = 0; blob_index < params.num_blobs; ++blob_index) {

		cout << "Blob " << blob_index << ":" << endl << endl;

		// Ignore if we have a static blob
		if(active_blob_array[blob_index]->get_motion_state() == FFEA_BLOB_IS_STATIC) {
			cout << "Blob " << blob_index << " is STATIC. No associated timesteps." << endl;
			continue;
		}

		// Define and reset all required variables for this crazy thing (maybe don't use stack memory??)
		num_nodes = active_blob_array[blob_index]->get_num_linear_nodes();
		num_rows = 3 * num_nodes;

		Eigen::SparseMatrix<scalar> K(num_rows, num_rows);
		Eigen::SparseMatrix<scalar> A(num_rows, num_rows);
		Eigen_MatrixX K_inv(num_rows, num_rows);
		Eigen_MatrixX I(num_rows, num_rows);
		Eigen_MatrixX tau_inv(num_rows, num_rows);

		// Build K
		cout << "\tCalculating the Global Viscosity Matrix, K...";
		if(active_blob_array[blob_index]->build_linear_node_viscosity_matrix(&K) == FFEA_ERROR) {
			cout << endl;
			FFEA_error_text();
			cout << "In function 'Blob::build_linear_node_viscosity_matrix' from blob " << blob_index << endl;
			return FFEA_ERROR;
		}
		cout << "done!" << endl;

		// Build A
		cout << "\tCalculating the Global Linearised Elasticity Matrix, A...";
		if(active_blob_array[blob_index]->build_linear_node_elasticity_matrix(&A) == FFEA_ERROR) {
			cout << endl;
			FFEA_error_text();
			cout << "In function 'Blob::build_linear_node_elasticity_matrix'" << blob_index << endl;
			return FFEA_ERROR;
		}
		cout << "done!" << endl;

		// Invert K (it's symmetric! Will not work if stokes_visc == 0)
		cout << "\tAttempting to invert K to form K_inv...";
		Eigen::SimplicialCholesky<Eigen::SparseMatrix<scalar>> Cholesky(K); // performs a Cholesky factorization of K
		I.setIdentity();
		K_inv = Cholesky.solve(I);
		if(Cholesky.info() == Eigen::Success) {
			cout << "done! Successful inversion of K!" << endl;
		} else if (Cholesky.info() == Eigen::NumericalIssue) {
			cout << "\nK cannot be inverted via Cholesky factorisation due to numerical issues. You possible don't have an external solvent set." << endl;
			return FFEA_OK;
		} else if (Cholesky.info() == Eigen::NoConvergence) {
			cout << "\nInversion iteration couldn't converge. K must be a crazy matrix. Possibly has zero eigenvalues?" << endl;
			return FFEA_OK;
		}

		// Apply to A
		cout << "\tCalculating inverse time constant matrix, tau_inv = K_inv * A...";
		tau_inv = K_inv * A;
		cout << "done!" << endl;

		// Diagonalise
		cout << "\tDiagonalising tau_inv..." << flush;
		Eigen::EigenSolver<Eigen_MatrixX> es(tau_inv);


		cout << "done!" << endl;

		// Output. Ignore the 6 translational modes (they are the slowest 6 modes)
		scalar smallest_val = fabs(es.eigenvalues()[6].real());
		scalar largest_val = fabs(es.eigenvalues()[0].real());

		cout << "\tThe time-constant of the slowest mode in Blob " << blob_index << ", tau_max = " << (1.0 / smallest_val) * mesoDimensions::time << "s" << endl;
		cout << "\tThe time-constant of the fastest mode in Blob " << blob_index << ", tau_min = " << (1.0 / largest_val) * mesoDimensions::time << "s" << endl << endl;

		// Global stuff
		if(1.0 / smallest_val > dt_max) {
			dt_max = 1.0 / smallest_val;
			dt_max_bin = blob_index;
		}

		if(1.0 / largest_val < dt_min) {
			dt_min = 1.0 / largest_val;
			dt_min_bin = blob_index;
		}

	}
	cout << "The time-constant of the slowest mode in all blobs, tau_max = " << dt_max * mesoDimensions::time << "s, from Blob " << dt_max_bin << endl;
	cout << "The time-constant of the fastest mode in all blobs, tau_min = " << dt_min * mesoDimensions::time << "s, from Blob " << dt_min_bin << endl << endl;
	cout << "Remember, the energies in your system will begin to be incorrect long before the dt = tau_min. I'd have dt << tau_min if I were you." << endl;
	return FFEA_OK;
}*/

/**
 * @brief Calculates an linear elastic model for a given blob.
 * @param[in] set<int> List of blobs to get Linear Elastic Model
 * @param[in] int num_modes The number of modes / eigenvalues to be calculated
 * @details By linearising the elasticity vector around the initial position,
 * this function performs matrix algebra using the Eigen libraries to diagonalise
 * elasticity matrix and output pseudo-trajectories based upon these eigenvectors
 * */
int World::lem(set<int> blob_indices, int num_modes)
{

    int i, j;
    int num_nodes, num_rows;
    set<int>::iterator it;

    vector<string> all;
    string traj_out_fname, base, ext, evals_out_fname, evecs_out_fname;

    // For all blobs in this set, calculate and output the elastic normal modes
    for (it = blob_indices.begin(); it != blob_indices.end(); it++)
    {

        // Get the index
        i = *it;

        cout << "\tBlob " << i << ":" << endl
             << endl;

        // Get an elasticity matrix
        num_nodes = active_blob_array[i]->get_num_linear_nodes();

        // Check this is allowed
        if (num_modes > 3 * num_nodes - 6)
        {
            cout << "\n\t\t" << num_modes << " unavailable for only " << num_nodes << " linear nodes. Deafulting to 3N-6 = " << 3 * num_nodes - 6 << " modes." << endl
                 << endl;
            num_modes = 3 * num_nodes - 6;
        }
        num_rows = num_nodes * 3;
        Eigen::SparseMatrix<scalar> A(num_rows, num_rows);

        cout << "\t\tCalculating the Linearised Elasticity Matrix, A...";
        if (active_blob_array[i]->build_linear_node_elasticity_matrix(&A) == FFEA_ERROR)
        {
            cout << endl;
            FFEA_error_text();
            cout << "In function 'Blob::build_linear_node_elasticity_matrix'" << i << endl;
            return FFEA_ERROR;
        }
        cout << "done!" << endl;

        // Diagonalise to find the elastic modes
        cout << "\t\tDiagonalising A...";
        Eigen::SelfAdjointEigenSolver<Eigen_MatrixX> es(A);
        cout << "done!" << endl;

        // This matrix 'should' contain 6 zero modes, and then num_rows - 6 actual floppy modes
        // The most important mode corresponds to the smallest non-zero eigenvalue

        // Most important mode will have motion ~ largest system size. Get a length...
        scalar dx = -1 * INFINITY;
        vector3 min, max;
        active_blob_array[i]->get_min_max(&min, &max);
        if (max.x - min.x > dx)
        {
            dx = max.x - min.x;
        }
        if (max.y - min.y > dx)
        {
            dx = max.y - min.y;
        }
        if (max.z - min.z > dx)
        {
            dx = max.z - min.z;
        }

        dx /= 5.0;

        // Sort evals into correct units (N/m)
        scalar unitscaler = mesoDimensions::force / mesoDimensions::length;

        // Make some trajectories (ignoring the first 6)
        cout << "\t\tMaking trajectories from eigenvectors..." << endl;

        // Get filename base
        ostringstream bi;
        bi << i;
        boost::split(all, params.trajectory_out_fname, boost::is_any_of("."));
        ext = "." + all.at(all.size() - 1);
        base = boost::erase_last_copy(params.trajectory_out_fname, ext);
        for (j = 6; j < 6 + num_modes; ++j)
        {
            cout << "\t\t\tEigenvector " << j << " / Mode " << j - 6 << "...";

            // Get a filename end
            ostringstream mi;
            mi << j - 6;
            traj_out_fname = base + "_FFEAlem_blob" + bi.str() + "mode" + mi.str() + ext;
            make_trajectory_from_eigenvector(traj_out_fname, i, j - 6, es.eigenvectors().col(j), dx);
            cout << "done!" << endl;
        }
        cout << "\t\tdone!" << endl;

        // Print out relevant eigenvalues and eigenvectors

        // Get a filename
        evals_out_fname = base + "_FFEAlem_blob" + bi.str() + ".evals";
        evecs_out_fname = base + "_FFEAlem_blob" + bi.str() + ".evecs";

        print_evecs_to_file(evecs_out_fname, es.eigenvectors(), num_rows, num_modes);
        print_evals_to_file(evals_out_fname, es.eigenvalues(), num_modes, unitscaler);
    }

    return FFEA_OK;
}

/**
 * @brief Calculates an dynamic mode model for a given blob.
 * @param[in] set<int> List of blobs to get Dynamic Mode Model
 * @param[in] int num_modes The number of modes / eigenvalues to be calculated
 * @details By linearising the elasticity vector around the initial position,
 * this function performs matrix algebra using the Eigen libraries to diagonalise
 * the coupling between the elasticity and viscosity matrices, and outputs
 * pseudo-trajectories based upon these eigenvectors
 * */

int World::dmm(set<int> blob_indices, int num_modes)
{

    int i, j, k;
    int num_nodes, num_rows;
    set<int>::iterator it;

    vector<string> all;
    string traj_out_fname, base, ext, evecs_out_fname, evals_out_fname;

    // For all blobs in this set, calculate and output the dynamic normal modes
    for (it = blob_indices.begin(); it != blob_indices.end(); it++)
    {

        // Get the index
        i = *it;

        cout << "\tBlob " << i << ":" << endl
             << endl;

        // Get a viscosity matrix
        num_nodes = active_blob_array[i]->get_num_linear_nodes();
        num_rows = num_nodes * 3;

        // Check this is allowed
        if (num_modes > 3 * num_nodes - 6)
        {
            cout << "\n\t\t" << num_modes << " unavailable for only " << num_nodes << " linear nodes. Deafulting to 3N-6 = " << 3 * num_nodes - 6 << " modes." << endl
                 << endl;
            num_modes = 3 * num_nodes - 6;
        }

        Eigen::SparseMatrix<scalar> K(num_rows, num_rows);

        cout << "\t\tCalculating the Global Linearised Viscosity Matrix, K...";
        if (active_blob_array[i]->build_linear_node_viscosity_matrix(&K) == FFEA_ERROR)
        {
            cout << endl;
            FFEA_error_text();
            cout << "In function 'Blob::build_linear_node_viscosity_matrix'" << i << endl;
            return FFEA_ERROR;
        }
        cout << "done!" << endl;

        // Diagonalise the thing
        cout << "\t\tDiagonalising K...";
        Eigen::SelfAdjointEigenSolver<Eigen_MatrixX> esK(K);
        cout << "done!" << endl;

        // Use this diagonalisation to define Q
        cout << "\t\tBuilding the matrix Q from the eigenvalues of K...";
        Eigen::SparseMatrix<scalar> Q(num_rows, num_rows);
        std::vector<Eigen::Triplet<scalar>> vals;
        for (j = 0; j < num_rows; ++j)
        {
            vals.push_back(Eigen::Triplet<scalar>(j, j, 1.0 / sqrt(fabs(esK.eigenvalues()[j]))));
        }
        Q.setFromTriplets(vals.begin(), vals.end());
        cout << "done!" << endl;

        // Get an elasticity matrix
        Eigen::SparseMatrix<scalar> A(num_rows, num_rows);
        cout << "\t\tCalculating the Global Linearised Elasticity Matrix, A...";
        if (active_blob_array[i]->build_linear_node_elasticity_matrix(&A) == FFEA_ERROR)
        {
            cout << endl;
            FFEA_error_text();
            cout << "In function 'Blob::build_linear_node_elasticity_matrix'" << i << endl;
            return FFEA_ERROR;
        }
        cout << "done!" << endl;

        // From A, build the transformation Ahat
        cout << "\t\tBuilding the Transformation Matrix Ahat...";
        Eigen_MatrixX Ahat(num_rows, num_rows);
        Ahat = Q.transpose() * esK.eigenvectors().transpose() * A * esK.eigenvectors() * Q;
        cout << "done" << endl;

        // Diagonalise to find the dynamic modes
        cout << "\t\tDiagonalising Ahat...";
        Eigen::SelfAdjointEigenSolver<Eigen_MatrixX> esAhat(Ahat);
        cout << "done!" << endl;
        cout << "Building the Dynamic Modes Matrix R...";
        Eigen_MatrixX R;
        R = esK.eigenvectors() * Q * esAhat.eigenvectors();

        // This matrix 'should' contain 6 zero modes, and then num_rows - 6 actual floppy modes
        // The most important mode corresponds to the smallest non-zero eigenvalue

        // Most important mode will have motion ~ largest system size. Get a length...
        scalar dx = -1 * INFINITY;
        vector3 min, max;
        active_blob_array[i]->get_min_max(&min, &max);
        if (max.x - min.x > dx)
        {
            dx = max.x - min.x;
        }
        if (max.y - min.y > dx)
        {
            dx = max.y - min.y;
        }
        if (max.z - min.z > dx)
        {
            dx = max.z - min.z;
        }

        dx /= 10.0;

        // Make some trajectories (ignoring the first 6)
        // Firstly, normalise the first num_modes eigenvectors
        scalar sum;
        for (j = 6; j < num_modes + 6; ++j)
        {
            sum = 0;
            for (k = 0; k < num_rows; ++k)
            {
                sum += R.col(j)[k] * R.col(j)[k];
            }
            R.col(j) *= 1.0 / sqrt(sum);
        }

        cout << "\t\tMaking trajectories from eigenvectors..." << endl;

        // Get filename base
        ostringstream bi, mi;
        bi << i;
        boost::split(all, params.trajectory_out_fname, boost::is_any_of("."));
        ext = "." + all.at(all.size() - 1);
        base = boost::erase_last_copy(params.trajectory_out_fname, ext);

        for (j = 6; j < 6 + num_modes; ++j)
        {

            cout << "\t\t\tEigenvector " << j << " / Mode " << j - 6 << "...";

            // Get a filename end
            ostringstream mi;
            mi << j - 6;
            traj_out_fname = base + "_FFEAdmm_blob" + bi.str() + "mode" + mi.str() + ext;

            make_trajectory_from_eigenvector(traj_out_fname, i, j - 6, R.col(j), dx);
            cout << "done!" << endl;
        }
        cout << "\t\tdone!" << endl;

        // Print out relevant eigenvalues and eigenvectors

        // Get a filename
        evals_out_fname = base + "_FFEAdmm_blob" + bi.str() + ".evals";
        evecs_out_fname = base + "_FFEAdmm_blob" + bi.str() + ".evecs";

        print_evecs_to_file(evecs_out_fname, R, num_rows, num_modes);
        print_evals_to_file(evals_out_fname, esAhat.eigenvalues(), num_modes, 1.0);
    }
    return FFEA_OK;
}

/**
 * @brief Calculates an elastic network model for a given blob.
 * @param[in] set<int> List of blobs to get Elastic Network Model
 * @param[in] int num_modes The number of modes / eigenvalues to be calculated
 * @details By linearising the elasticity vector around the initial position,
 * this function performs matrix algebra using the Eigen libraries to diagonalise
 * the coupling between the elasticity and rotne-praga viscosity matrices, and outputs
 * pseudo-trajectories based upon these eigenvectors
 * */

int World::dmm_rp(set<int> blob_indices, int num_modes)
{

    int i, j, k, l;
    int num_nodes, num_rows;
    set<int>::iterator it;

    vector<string> all;
    string traj_out_fname, base, ext, evecs_out_fname, evals_out_fname;

    // For all blobs in this set, calculate and output the dynamic normal modes
    for (it = blob_indices.begin(); it != blob_indices.end(); it++)
    {

        // Get the index
        i = *it;

        cout << "\tBlob " << i << ":" << endl
             << endl;

        // Explicitly calculate a diffusion matrix
        num_nodes = active_blob_array[i]->get_num_linear_nodes();
        num_rows = num_nodes * 3;

        // Check this is allowed
        if (num_modes > 3 * num_nodes - 6)
        {
            cout << "\n\t\t" << num_modes << " unavailable for only " << num_nodes << " linear nodes. Deafulting to 3N-6 = " << 3 * num_nodes - 6 << " modes." << endl
                 << endl;
        }
        num_modes = 3 * num_nodes - 6;

        Eigen_MatrixX D(num_rows, num_rows);

        cout << "\t\tCalculating the Rotne-Prager diffusion matrix, D..." << flush;
        if (active_blob_array[i]->build_linear_node_rp_diffusion_matrix(&D) == FFEA_ERROR)
        {
            cout << endl;
            FFEA_error_text();
            cout << "In function 'Blob::build_rp_diffusion_matrix'" << i << endl;
            return FFEA_ERROR;
        }
        cout << "done!" << endl;

        // Get an elasticity matrix
        Eigen::SparseMatrix<scalar> A(num_rows, num_rows);
        cout << "\t\tCalculating the Global Linearised Elasticity Matrix, A..." << flush;
        if (active_blob_array[i]->build_linear_node_elasticity_matrix(&A) == FFEA_ERROR)
        {
            cout << endl;
            FFEA_error_text();
            cout << "In function 'Blob::build_linear_node_elasticity_matrix'" << i << endl;
            return FFEA_ERROR;
        }
        cout << "done!" << endl;

        // Diagonalise DA to find the dynamic modes
        cout << "\t\tCalculating D*A..." << flush;
        Eigen_MatrixX F;
        F = D * A;
        cout << "done!" << endl;
        cout << "\t\tDiagonalising DA..." << flush;
        Eigen::EigenSolver<Eigen_MatrixX> esF(F);
        cout << "done!" << endl;

        // Order the eigenvalues
        Eigen_MatrixX Rvecs(num_rows, num_modes + 6);
        Eigen_VectorX Rvals(num_modes + 6);
        Eigen::VectorXi Rvals_indices(num_modes + 6);
        Rvals.setZero();
        Rvecs.setZero();
        scalar min_eig, max_eig, aneig;
        int max_eig_index, fail;

        // Get the 6 zero modes and the num_modes required modes in order
        for (j = 0; j < num_modes + 6; ++j)
        {

            // Set limits
            max_eig = INFINITY;
            if (j == 0)
            {
                min_eig = -1 * INFINITY;
            }
            else
            {
                min_eig = Rvals(j - 1);
            }

            // Sweep the eigenvalue range
            for (k = 0; k < num_rows; ++k)
            {

                // Ignore imaginary stuff (should be none anyway, but you know how numerical errors are...)
                if (esF.eigenvalues()[k].imag() != 0)
                {
                    continue;
                }
                aneig = fabs(esF.eigenvalues()[k].real());

                // If within limits
                if (aneig <= max_eig && aneig >= min_eig)
                {

                    // Check if already exists in vector (probably quicker to use a set but I've started now!)
                    fail = 0;
                    for (l = 0; l < j; ++l)
                    {
                        if (Rvals_indices(l) == k)
                        {
                            fail = 1;
                            break;
                        }
                    }

                    // Change limits
                    if (fail == 0)
                    {
                        max_eig = aneig;
                        max_eig_index = k;
                    }
                }
            }
            Rvals_indices(j) = max_eig_index;
            Rvals(j) = esF.eigenvalues()[max_eig_index].real();
            Rvecs.col(j) = esF.eigenvectors().col(max_eig_index).real();
        }

        // This matrix 'should' contain 6 zero modes, and then num_rows - 6 actual floppy modes
        // The most important mode corresponds to the smallest non-zero eigenvalue

        // Most important mode will have motion ~ largest system size. Get a length...
        scalar dx = -1 * INFINITY;
        vector3 min, max;
        active_blob_array[i]->get_min_max(&min, &max);
        if (max.x - min.x > dx)
        {
            dx = max.x - min.x;
        }
        if (max.y - min.y > dx)
        {
            dx = max.y - min.y;
        }
        if (max.z - min.z > dx)
        {
            dx = max.z - min.z;
        }

        dx /= 20.0;

        // Make some trajectories (ignoring the first 6)
        // Firstly, normalise the first num_modes eigenvectors
        scalar sum;
        for (j = 6; j < num_modes + 6; ++j)
        {
            sum = 0;
            for (k = 0; k < num_rows; ++k)
            {
                sum += Rvecs.col(j)[k] * Rvecs.col(j)[k];
            }
            Rvecs.col(j) *= 1.0 / sqrt(sum);
        }

        cout << "\t\tMaking trajectories from eigenvectors..." << endl;

        // Get filename base
        ostringstream bi, mi;
        bi << i;
        boost::split(all, params.trajectory_out_fname, boost::is_any_of("."));
        ext = "." + all.at(all.size() - 1);
        base = boost::erase_last_copy(params.trajectory_out_fname, ext);

        for (j = 6; j < 6 + num_modes; ++j)
        {

            cout << "\t\t\tEigenvector " << j << " / Mode " << j - 6 << "...";

            // Get a filename end
            ostringstream mi;
            mi << j - 6;
            traj_out_fname = base + "_FFEArpdmm_blob" + bi.str() + "mode" + mi.str() + ext;

            make_trajectory_from_eigenvector(traj_out_fname, i, j - 6, Rvecs.col(j), dx);
            cout << "done!" << endl;
        }
        cout << "\t\tdone!" << endl;

        // Print out relevant eigenvalues and eigenvectors

        // Get a filename
        evals_out_fname = base + "_FFEArpdmm_blob" + bi.str() + ".evals";
        evecs_out_fname = base + "_FFEArpdmm_blob" + bi.str() + ".evecs";

        print_evecs_to_file(evecs_out_fname, Rvecs, num_rows, num_modes);
        print_evals_to_file(evals_out_fname, Rvals, num_modes, 1.0);
    }
    return FFEA_OK;
}

/**
 * Update entire World for num_steps time steps
 * */
int World::run()
{

    int es_count = params.es_update;
    scalar wtime = omp_get_wtime();
#ifdef BENCHMARK
    scalar wtime0, wtime1, wtime2, wtime3, wtime4, time0 = 0, time1 = 0, time2 = 0, time3 = 0, time4 = 0;
#endif

    for (long long step = step_initial; step < params.num_steps + 1; step++)
    {
#ifdef BENCHMARK
        wtime0 = omp_get_wtime();
#endif

        // check if we are in one of those time-steps
        //     where we need to calculate electrostatics,
        //                     update the neighbour list,
        //                     and calculate normals and centroids for the faces.
        bool es_update = false;
        if (params.calc_ssint == 1 || params.calc_steric == 1 || params.calc_preComp == 1 || params.sticky_wall_xz == 1 || params.calc_es == 1)
        {
            if (es_count == params.es_update)
            {
                es_count = 1;
                es_update = true;
#ifdef FFEA_PARALLEL_FUTURE
                // we will have to wait for both threads to finish!
                if (updatingVdWLL() == true)
                {
                    cout << "Waiting for thread VdW LL" << endl;
                    if (catch_thread_updatingVdWLL(step, wtime, 1))
                        return FFEA_ERROR;
                }
                if (updatingPCLL() == true)
                {
                    cout << "Waiting for thread PC LL" << endl;
                    if (catch_thread_updatingPCLL(step, wtime, 1))
                        return FFEA_ERROR;
                }
#endif
            }
            else
                es_count++;
        } // es_update will turn to false at the begining of next timestep

        // Zero the force across all blobs
#ifdef USE_OPENMP
#pragma omp parallel for default(none) shared(es_update, step) schedule(static)
#endif
        for (int i = 0; i < params.num_blobs; i++)
        {
            active_blob_array[i]->zero_force();
            /* if ((step+1) % params.check == 0) { // we only do measurements if we need so.
                // if (params.calc_vdw == 1) active_blob_array[i]->zero_vdw_bb_measurement_data(); // DEPRECATED
                // if (params.sticky_wall_xz == 1) active_blob_array[i]->zero_vdw_xz_measurement_data(); // DEPRECATED
            } */

            // If blob centre of mass moves outside simulation box, apply PBC to it
            vector3 com;
            active_blob_array[i]->calc_and_store_centroid(com);

            scalar dx = 0, dy = 0, dz = 0;
            int check_move = 0;

            if (com.x < 0)
            {
                if (params.wall_x_1 == WALL_TYPE_PBC)
                {
                    dx += box_dim.x;
                    active_blob_array[i]->pbc_count[0] -= 1;
                    check_move = 1;
                }
            }
            else if (com.x > box_dim.x)
            {
                if (params.wall_x_2 == WALL_TYPE_PBC)
                {
                    dx -= box_dim.x;
                    active_blob_array[i]->pbc_count[0] += 1;
                    check_move = 1;
                }
            }
            if (com.y < 0)
            {
                if (params.wall_y_1 == WALL_TYPE_PBC)
                {
                    dy += box_dim.y;
                    active_blob_array[i]->pbc_count[1] -= 1;
                    check_move = 1;
                }
            }
            else if (com.y > box_dim.y)
            {
                if (params.wall_y_2 == WALL_TYPE_PBC)
                {
                    dy -= box_dim.y;
                    active_blob_array[i]->pbc_count[1] += 1;
                    check_move = 1;
                }
            }
            if (com.z < 0)
            {
                if (params.wall_z_1 == WALL_TYPE_PBC)
                {
                    dz += box_dim.z;
                    active_blob_array[i]->pbc_count[2] -= 1;
                    check_move = 1;
                }
            }
            else if (com.z > box_dim.z)
            {
                if (params.wall_z_2 == WALL_TYPE_PBC)
                {
                    dz -= box_dim.z;
                    active_blob_array[i]->pbc_count[2] += 1;
                    check_move = 1;
                }
            }
            // So if it has to move, do so and update its centroid.
            if (check_move == 1)
            {
                active_blob_array[i]->move(dx, dy, dz);
                active_blob_array[i]->calc_and_store_centroid(com);
            }

            // If Blob is near a hard wall, prevent it from moving further into it
            active_blob_array[i]->enforce_box_boundaries(&box_dim);

            // Set node forces to zero
            active_blob_array[i]->set_forces_to_zero();

            if (es_update)
                active_blob_array[i]->calc_centroids_and_normals_of_all_faces();
        }

        if (es_update)
        {
            // REFRESH LINKED LISTS:
            // Refresh the VdW-LinkedList
#ifdef FFEA_PARALLEL_FUTURE
            // Thread out to update the LinkedLists,
            //   after calculating the centroids of the faces.
            // Catching up the thread should be done through catch_thread_updatingVdWLL,
            //   which will to lookup.safely_swap_layers().
            if (updatingVdWLL() == false)
            {
                thread_updatingVdWLL = std::async(std::launch::async, &World::prebuild_nearest_neighbour_lookup_wrapper, this, params.ssint_cutoff);
            }
            else // die with dignity
#else
            if (lookup.build_nearest_neighbour_lookup(params.ssint_cutoff) == FFEA_ERROR)
#endif
            {
                return die_with_dignity(step, wtime);
            }

            // Refresh the PreComp-LinkedList
#ifdef FFEA_PARALLEL_FUTURE
            // Just as in VdW, thread out to update the PC LinkedLists
            if (updatingPCLL() == false)
            {
                thread_updatingPCLL = std::async(std::launch::async, &PreComp_solver::prebuild_pc_nearest_neighbour_lookup, &pc_solver);
            }
            else // die with dignity
#else
            if (pc_solver.build_pc_nearest_neighbour_lookup() == FFEA_ERROR)
#endif
            {
                return die_with_dignity(step, wtime);
            }
            // FINSHED REFRESHING LINKED LISTS.

            // Finally do calc_es, which is done only from time to time...
            if (params.calc_es == 1)
                do_es();
        }

        // Apply springs directly to nodes
        apply_springs();

        // Rod steric neighbours
        if (params.calc_steric_rod == 1)
        {
            for (int i = 0; i < params.num_rods; i++)
            {
                for (int j = i + 1; j < params.num_rods; j++)
                {
                    update_rod_steric_nbr_lists(rod_array[i], rod_array[j]);
                }
            }
            if (rod::dbg_print)
            {
                std::cout << "Generated rod steric neighbour lists" << std::endl;
            }
        }
        else if (rod::dbg_print)
        {
            std::cout << "Rod-rod steric interactions disabled." << std::endl;
        }

        // Rod VDW neighbours
        if (params.calc_vdw_rod == 1)
        {
            for (int i = 0; i < params.num_rods; i++)
            {
                for (int j = i + 1; j < params.num_rods; j++)
                {
                    update_rod_vdw_nbr_lists(rod_array[i], rod_array[j], &rod_lj_matrix);
                }
            }
            if (rod::dbg_print)
            {
                std::cout << "Generated rod vdw neighbour lists" << std::endl;
            }
        }
        else if (rod::dbg_print)
        {
            std::cout << "Rod-rod vdw interactions disabled." << std::endl;
        }


        // Do rods
        for (int i = 0; i < params.num_rods; i++)
        {
            rod_array[i]->do_timestep(rng);
        }

        // Rod-blob interface
        for (int i = 0; i < params.num_interfaces; i++)
        {
            rod_blob_interface_array[i]->do_connection_timestep();
        }

        if (params.pbc_rod == 1)
        {
            for (int i = 0; i < params.num_rods; i++)
                rod_pbc_wrap(rod_array[i], {(float)box_dim.x, (float)box_dim.y, (float)box_dim.z});
        }

#ifdef FFEA_PARALLEL_FUTURE
        // Get the thread updating the VdW-LinkedLists if it has finished.
        if (updatingPCLL_ready_to_swap() == true)
        {
            if (catch_thread_updatingPCLL(step, wtime, 3))
                return die_with_dignity(step, wtime);
        }
#endif

#ifdef BENCHMARK
        wtime1 = omp_get_wtime();
        time0 += wtime1 - wtime0;
#endif

        // Calculate the correlation matrix if Blob-Blob forces need PBC:
        if (params.force_pbc == 1)
            calc_blob_corr_matrix(params.num_blobs, blob_corr);

        // Calculate the PreComp forces:
        if (params.calc_preComp == 1)
        {
            // pc_solver.solve(blob_corr); // keep that one as reference implementation.
            // pc_solver.solve_using_neighbours();
            pc_solver.solve_using_neighbours_non_critical(blob_corr); // blob_corr == NULL if force_pbc = 0
        }

#ifdef BENCHMARK
        wtime2 = omp_get_wtime();
        time1 += wtime2 - wtime1;
#endif

#ifdef FFEA_PARALLEL_FUTURE
        // Get the thread updating the VdW-LinkedLists if it has finished.
        // #pragma omp master // Then a single thread does the catching and swapping
        if (updatingVdWLL_ready_to_swap() == true)
        {
            if (catch_thread_updatingVdWLL(step, wtime, 3))
                return die_with_dignity(step, wtime);
        }
        // #pragma omp barrier // the barrier holds people off, before catching the thread
#endif

        // Calculate the VdW forces:
        if (params.calc_ssint == 1 || params.calc_steric == 1)
        {
            vdw_solver->solve(blob_corr); // blob_corr == NULL if force_pbc = 0.
        }
        if (params.sticky_wall_xz == 1)
            vdw_solver->solve_sticky_wall(params.ssint_cutoff);

            //checks whether force periodic boundary conditions specified, calculates periodic array correction to array through vdw_solver as overload

#ifdef BENCHMARK
        wtime3 = omp_get_wtime();
        time2 += wtime3 - wtime2;
#endif

        // Output to checkpoint files the random numbers responsible for the current system before they advance in the update_internal_forces bit
        if ((step) % params.check == 0)
        {
            print_checkpoints();
        }

        // Update Blobs, while tracking possible errors.
        int fatal_errors = 0;
#ifdef FFEA_PARALLEL_PER_BLOB
#pragma omp parallel for default(none) reduction(+ \
                                                 : fatal_errors) schedule(static)
#endif
        for (int i = 0; i < params.num_blobs; i++)
        {
            if (active_blob_array[i]->update_internal_forces() == FFEA_ERROR)
                fatal_errors++;
        }

        if (fatal_errors > 0)
            return die_with_dignity(step, wtime);

        // Now, for this step, all forces and positions are correct, as are the kinetic states
        if ((step) % params.check == 0)
        {
            print_trajectory_and_measurement_files(step, wtime);
            print_kinetic_files(step);
        }

        // Finally, update the positions
#ifdef FFEA_PARALLEL_PER_BLOB

#pragma omp parallel for default(none) schedule(static)
#endif
        for (int i = 0; i < params.num_blobs; i++)
        {
            active_blob_array[i]->update_positions();
        }

#ifdef BENCHMARK
        wtime4 = omp_get_wtime();
        time3 += wtime4 - wtime3;
#endif

        /* Kinetic Part of each step */
        // This part consists of a discrete change, and so must occur before a force calculation cycle to be consistent with measurement data
        // This means is must happen either at the very end of a timstep, or at the very beginning
        if (params.calc_kinetics == 1 && step % params.kinetics_update == 0)
        {

            // Calculate the kinetic switching probablilites. These are scaled from the base rates provided
            if (calculate_kinetic_rates() != FFEA_OK)
            {
                FFEA_ERROR_MESSG("'calculate_kinetic_rates()' failed.\n")
            }

            // Now we can treat each blob separately
            int target;
            for (int i = 0; i < params.num_blobs; ++i)
            {

                // Find out what the state change should be based on the rates
                target = active_blob_array[i]->get_state_index();
                if (choose_new_kinetic_state(i, &target) != FFEA_OK)
                {
                    FFEA_ERROR_MESSG("'calculate_kinetic_rates()' failed.\n")
                }

                if (change_kinetic_state(i, target) != FFEA_OK)
                {
                    FFEA_ERROR_MESSG("'change_kinetic_state()' failed.\n")
                }
            }
        }

#ifdef BENCHMARK
        time4 += omp_get_wtime() - wtime4;
#endif
    }

#ifdef FFEA_PARALLEL_FUTURE
    // Wait until the last step has correctly been written:
    thread_writingTraj.get();
#endif

    // Total mpi timing, compare with openmp timing
#ifdef BENCHMARK
    cout << "total steps:" << params.num_steps << endl;
    cout << "benchmarking--------cleaning & threading out\t" << time0 << "seconds" << endl;
    cout << "benchmarking--------strictly preComp for \t" << time1 << "seconds" << endl;
    cout << "benchmarking--------ssint & threading back in\t" << time2 << "seconds" << endl;
    cout << "benchmarking--------update blobs for \t" << time3 << "seconds" << endl;
    cout << "benchmarking--------writing and measuring \t" << time4 << "seconds" << endl;
    cout << "benchmarking--------Total time in World::run():" << omp_get_wtime() - wtime << "seconds" << endl;

    /* cout<< "benchmarking--------PreComp decomposition: " << endl;
    cout<< "            --------Cleaning: " << pc_solver.timep0 << " seconds" << endl;
    cout<< "            --------Positions: " << pc_solver.timep1 << " seconds" << endl;
    cout<< "            --------Doubleloop: " << pc_solver.timep2 << " seconds" << endl;
    cout<< "            ------------------part1 : " << pc_solver.timep3 << " seconds" << endl;
    cout<< "            ------------------part2 : " << pc_solver.timep4 << " seconds" << endl; */
#endif

    printf("\n\nTime taken: %2f seconds\n", (omp_get_wtime() - wtime));

    return FFEA_OK;
}

/**
 * @brief Changes the active kinetic state for a given blob.
 * @param[in] int blob_index Index of the blob to be changed
 * @param[in] int target_state State the blob should be changed to
 * @details This function changes the kinetic state of the blob.
 * For a conformational change, the current structure is mapped to
 * the target and the active_blob_array[i] pointer is updated.
 * For a binding / unbinding event, the binding sites are activated.
 * */

int World::change_kinetic_state(int blob_index, int target_state)
{

    // Do we even need to change anything?
    int current_state = active_blob_array[blob_index]->get_state_index();
    if (current_state == target_state)
    {
        return FFEA_OK;
    }

    int current_conformation = active_blob_array[blob_index]->get_conformation_index();
    int target_conformation = kinetic_state[blob_index][target_state].get_conformation_index();

    // If we do, how do we change?
    if (kinetic_state[blob_index][current_state].get_conformation_index() != kinetic_state[blob_index][target_state].get_conformation_index())
    {

        // Conformation change!
        int inversionCheck;

        // Get current nodes
        vector3 **current_nodes = active_blob_array[blob_index]->get_actual_node_positions();

        // Get target nodes
        vector3 **target_nodes = blob_array[blob_index][target_conformation].get_actual_node_positions();

        // Apply map
        kinetic_map[blob_index][current_conformation][target_conformation].block_apply(current_nodes, target_nodes);

        // Check inversion
        inversionCheck = blob_array[blob_index][target_conformation].check_inversion();

        if (inversionCheck == FFEA_ERROR)
        {
            printf("Conformational change rejected.\n");

            // Move target conformation back to infinity
            blob_array[blob_index][target_conformation].position(blob_array[blob_index][target_conformation].get_RandU01() * 1e10, blob_array[blob_index][target_conformation].get_RandU01() * 1e10, blob_array[blob_index][target_conformation].get_RandU01() * 1e10);

            return FFEA_OK;
        }
        else
        {

            // Change active conformation and activate all faces
            active_blob_array[blob_index] = &blob_array[blob_index][target_conformation];
            active_blob_array[blob_index]->kinetically_set_faces(true);

            // Move the old one to random space so as not to interfere with calculations, and deactivate all faces
            blob_array[blob_index][current_conformation].position(blob_array[blob_index][current_conformation].get_RandU01() * 1e10, blob_array[blob_index][current_conformation].get_RandU01() * 1e10, blob_array[blob_index][current_conformation].get_RandU01() * 1e10);
            blob_array[blob_index][current_conformation].kinetically_set_faces(false);

            // Reactivate springs
            activate_springs();
        }
    }
    else if (!kinetic_state[blob_index][current_state].is_bound() && kinetic_state[blob_index][target_state].is_bound())
    {

        // Binding event! Add nodes to pinned node list, or add springs, and reset the solver
        active_blob_array[blob_index]->pin_binding_site(kinetic_state[blob_index][target_state].get_base_site()->get_nodes());
        active_blob_array[blob_index]->reset_solver();
    }
    else if (kinetic_state[blob_index][current_state].is_bound() && !kinetic_state[blob_index][target_state].is_bound())
    {

        // Unbinding event! Remove nodes to pinned node list and reset the solver
        active_blob_array[blob_index]->unpin_binding_site(kinetic_state[blob_index][current_state].get_base_site()->get_nodes());
        active_blob_array[blob_index]->reset_solver();
    }
    else
    {

        // Identity event. Nothing happens
    }

    // Change all indices
    active_blob_array[blob_index]->set_previous_state_index(current_state);
    active_blob_array[blob_index]->set_state_index(target_state);
    active_blob_array[blob_index]->set_previous_conformation_index(current_conformation);

    return FFEA_OK;
}

/**
 * @brief Parses <blobs>, <springs>, <rods> and <precomp>.
 * @param[in] vector<string> script_vector, which is essentially the FFEA input file,
 *            line by line, as it comes out of FFEA_input_reader::file_to_lines
 */
int World::read_and_build_system(vector<string> script_vector)
{

    // READ and parse more
    // Reading variables
    systemreader = new FFEA_input_reader();
    int i, j;
    string lrvalue[2]; //, maplvalue[2];
    vector<string> blob_vector, interactions_vector, conformation_vector, kinetics_vector, map_vector, spring_vector;
    vector<string>::iterator it;

    vector<string> nodes, topology, surface, material, stokes, ssint, binding, pin, beads;
    string map_fname;
    int map_indices[2];
    int set_motion_state = 0, set_nodes = 0, set_top = 0, set_surf = 0, set_mat = 0, set_stokes = 0, set_ssint = 0, set_binding = 0, set_pin = 0, set_solver = 0, set_preComp = 0, set_scale = 0, set_states = 0, set_rates = 0, calc_compress = 0;
    scalar scale = 1, compress = 1;
    int solver = FFEA_NOMASS_CG_SOLVER;
    vector<int> motion_state;

    vector<Blob_conf> blob_conf;

    // Read the INTERACTIONS block:
    // Get interactions vector first, for later use
    if ((params.calc_preComp == 1) || (params.calc_springs == 1) || (params.calc_ctforces == 1))
    {
        if (systemreader->extract_block("interactions", 0, script_vector, &interactions_vector) == FFEA_ERROR)
            return FFEA_ERROR;
    }

    // Get precomputed data first
    pc_params.dist_to_m = 1;
    pc_params.E_to_J = 1;
    if (params.calc_preComp == 1)
    {
        vector<string> precomp_vector;
        if (systemreader->extract_block("precomp", 0, interactions_vector, &precomp_vector) == FFEA_ERROR)
        {
            return FFEA_ERROR;
        }

        for (i = 0; i < precomp_vector.size(); i++)
        {
            systemreader->parse_tag(precomp_vector[i], lrvalue);
            if (lrvalue[0] == "types")
            {
                lrvalue[1] = boost::erase_last_copy(boost::erase_first_copy(lrvalue[1], "("), ")");
                boost::trim(lrvalue[1]);
                if (lrvalue[1].compare("") == 0)
                {
                    FFEA_ERROR_MESSG("Invalid value for 'types' in <precomp> section\n");
                    return FFEA_ERROR;
                }
                systemreader->split_string(lrvalue[1], pc_params.types, ",");
            }
            else if (lrvalue[0] == "inputData")
            {
                pc_params.inputData = stoi(lrvalue[1]);
            }
            else if (lrvalue[0] == "folder")
            {
                fs::path auxpath = params.FFEA_script_path / lrvalue[1];
                pc_params.folder = auxpath.string(); //   lrvalue[1];
            }
            else if (lrvalue[0] == "dist_to_m")
            {
                pc_params.dist_to_m = stod(lrvalue[1]);
            }
            else if (lrvalue[0] == "E_to_J")
            {
                pc_params.E_to_J = stod(lrvalue[1]);
            }
        }
    }

    if (params.calc_ctforces == 1)
    {
        vector<string> ctforces_vector;
        if (systemreader->extract_block("ctforces", 0, interactions_vector, &ctforces_vector) == FFEA_ERROR)
        {
            return FFEA_ERROR;
        }
        if (ctforces_vector.size() != 1)
        {
            FFEA_ERROR_MESSG("only a single field is allowed in <ctforces> block\n");
            return FFEA_ERROR;
        }
        systemreader->parse_tag(ctforces_vector[0], lrvalue);
        if (lrvalue[0] == "ctforces_fname")
        {
            fs::path auxpath = params.FFEA_script_path / lrvalue[1];
            params.ctforces_fname = auxpath.string();
        }
        else
        {
            FFEA_ERROR_MESSG("ctforces_fname not provided in <ctforces> block, while calc_ctforces set to 1\n");
            return FFEA_ERROR;
        }
        // try to open this file:
        if (!fs::exists(params.ctforces_fname))
        {
            FFEA_ERROR_MESSG("ctforces_fname: %s could not be open\n", params.ctforces_fname.c_str());
        }
    }

    if (params.calc_springs == 1)
    {
        if (systemreader->extract_block("springs", 0, interactions_vector, &spring_vector) == FFEA_ERROR)
            return FFEA_ERROR;

        if (spring_vector.size() > 1)
        {
            FFEA_error_text();
            cout << "'Spring' block should only have 1 file." << endl;
            return FFEA_ERROR;
        }
        else if (spring_vector.size() == 1)
        {
            systemreader->parse_tag(spring_vector.at(0), lrvalue);
            fs::path auxpath = params.FFEA_script_path / lrvalue[1];
            params.springs_fname = auxpath.string();
        }
    }
    params.write_to_file(stdout, pc_params);

    // DONE
    // Create some blobs based on params
    cout << "\tCreating blob array..." << endl;
    blob_array = new Blob *[params.num_blobs];
    active_blob_array = new Blob *[params.num_blobs];

#ifdef FFEA_PARALLEL_PER_BLOB
#pragma omp parallel for default(none) schedule(static) // shared(blob_array, active_blob_array)
#endif
    for (int i = 0; i < params.num_blobs; ++i)
    {
        blob_array[i] = new Blob[params.num_conformations[i]];
        active_blob_array[i] = &blob_array[i][0];
    }

    //creates blob correction array if specified in input file
    if (params.force_pbc == 1)
        blob_corr = new scalar[params.num_blobs * params.num_blobs * 3];

    // Read in each blob one at a time
    for (i = 0; i < params.num_blobs; ++i)
    {

        // Initialise some blob_conf stuff:
        blob_conf.push_back(Blob_conf());
        blob_conf[i].set_centroid = 0;
        blob_conf[i].set_velocity = 0;
        blob_conf[i].set_rotation = 0;

        // Get blob data
        if (systemreader->extract_block("blob", i, script_vector, &blob_vector) == FFEA_ERROR)
            return FFEA_ERROR;

        // Read all conformations
        bool enforce_conf_blocks = false;
        bool read_blob_as_conf = false;
        if (params.num_conformations[i] > 1)
            enforce_conf_blocks = true;
        for (j = 0; j < params.num_conformations[i]; ++j)
        {

            // Get conformation data
            int read_conf_err = systemreader->extract_block("conformation", j, blob_vector, &conformation_vector, enforce_conf_blocks);
            if (read_conf_err == FFEA_ERROR)
            {
                if (enforce_conf_blocks == true)
                    return FFEA_ERROR;
                else
                {
                    read_blob_as_conf = true;
                    conformation_vector = blob_vector;
                }
            }

            // Error check
            if (conformation_vector.size() == 0)
            {
                FFEA_error_text();
                cout << " In 'Blob' block " << i << ", expected at least a single 'conformation' block." << endl;
                return FFEA_ERROR;
            }

            // Parse conformation data
            for (it = conformation_vector.begin(); it != conformation_vector.end(); ++it)
            {
                systemreader->parse_tag(*it, lrvalue);

                // Assign if possible
                if (lrvalue[0] == "motion_state")
                {

                    if (lrvalue[1] == "DYNAMIC")
                    {
                        motion_state.push_back(FFEA_BLOB_IS_DYNAMIC);
                    }
                    else if (lrvalue[1] == "STATIC")
                    {
                        motion_state.push_back(FFEA_BLOB_IS_STATIC);
                    }
                    else if (lrvalue[1] == "FROZEN")
                    {
                        motion_state.push_back(FFEA_BLOB_IS_FROZEN);
                    }
                    set_motion_state = 1;
                }
                else if (lrvalue[0] == "nodes")
                {
                    fs::path auxpath = params.FFEA_script_path / lrvalue[1];
                    nodes.push_back(auxpath.string());
                    set_nodes = 1;
                }
                else if (lrvalue[0] == "topology")
                {
                    fs::path auxpath = params.FFEA_script_path / lrvalue[1];
                    topology.push_back(auxpath.string());
                    set_top = 1;
                }
                else if (lrvalue[0] == "surface")
                {
                    fs::path auxpath = params.FFEA_script_path / lrvalue[1];
                    surface.push_back(auxpath.string());
                    set_surf = 1;
                }
                else if (lrvalue[0] == "material")
                {
                    fs::path auxpath = params.FFEA_script_path / lrvalue[1];
                    material.push_back(auxpath.string());
                    set_mat = 1;
                }
                else if (lrvalue[0] == "stokes")
                {
                    fs::path auxpath = params.FFEA_script_path / lrvalue[1];
                    stokes.push_back(auxpath.string());
                    set_stokes = 1;
                }
                else if (lrvalue[0] == "ssint" || lrvalue[0] == "vdw")
                {
                    fs::path auxpath = params.FFEA_script_path / lrvalue[1];
                    ssint.push_back(auxpath.string());
                    set_ssint = 1;
                }
                else if (lrvalue[0] == "binding_sites")
                {
                    if (params.calc_kinetics == 1)
                    {
                        fs::path auxpath = params.FFEA_script_path / lrvalue[1];
                        binding.push_back(auxpath.string());
                        set_binding = 1;
                    }
                }
                else if (lrvalue[0] == "pin")
                {
                    fs::path auxpath = params.FFEA_script_path / lrvalue[1];
                    pin.push_back(auxpath.string());
                    set_pin = 1;
                }
                else if (lrvalue[0] == "beads")
                {
                    fs::path auxpath = params.FFEA_script_path / lrvalue[1];
                    beads.push_back(auxpath.string());
                    set_preComp = 1;
                }
                else
                {
                    bool forgive_unrecognised = false;
                    if ((read_blob_as_conf == true) && ((lrvalue[0] == "centroid") || (lrvalue[0] == "rotation") || (lrvalue[0] == "solver") || (lrvalue[0] == "scale") || (lrvalue[0] == "translate") || (lrvalue[0] == "velocity") || (lrvalue[0] == "calc_compress") || (lrvalue[0] == "compress")))
                        forgive_unrecognised = true;
                    if (forgive_unrecognised == false)
                    {
                        FFEA_error_text();
                        cout << "Unrecognised conformation lvalue: '" << lrvalue[0] << "'" << endl;
                        return FFEA_ERROR;
                    }
                }
            }

            // Error check
            if (set_motion_state == 0)
            {
                FFEA_error_text();
                cout << "Blob " << i << ", conformation " << j << " must have a motion state set.\nAccepted states: DYNAMIC, STATIC, FROZEN." << endl;
                return FFEA_ERROR;
            }
            else
            {
                if (set_nodes == 0 || set_surf == 0 || set_ssint == 0)
                {
                    FFEA_error_text();
                    cout << "In blob " << i << ", conformation " << j << ":\nFor any blob conformation, 'nodes', 'surface' and 'ssint' must be set." << endl;
                    return FFEA_ERROR;
                }

                if (motion_state.back() == FFEA_BLOB_IS_DYNAMIC)
                {
                    if (set_top == 0 || set_mat == 0 || set_stokes == 0)
                    {
                        FFEA_error_text();
                        cout << "In blob " << i << ", conformation " << j << ":\nFor a DYNAMIC blob conformation, 'topology', 'material', and 'stokes' must be set." << endl;
                        return FFEA_ERROR;
                    }
                }
                else
                {
                    topology.push_back("");
                    material.push_back("");
                    stokes.push_back("");
                    pin.push_back("");
                    set_top = 1;
                    set_mat = 1;
                    set_stokes = 1;
                    set_pin = 1;
                }

                // Optional stuff
                if (set_preComp == 0)
                {
                    beads.push_back("");
                    set_preComp = 1;
                }

                if (set_binding == 0)
                {
                    binding.push_back("");
                    set_binding = 1;
                }

                if (set_pin == 0)
                {
                    pin.push_back("");
                    set_pin = 1;
                }
            }

            // Clear conformation vector and set values for next round
            set_nodes = 0;
            set_top = 0;
            set_surf = 0;
            set_mat = 0;
            set_stokes = 0;
            set_ssint = 0;
            set_binding = 0;
            set_pin = 0;
            set_preComp = 0;
            conformation_vector.clear();
        }

        // Read kinetic info (if necessary)
        if (params.calc_kinetics == 1 && params.num_states[i] != 1)
        {

            // Get kinetic data
            if (systemreader->extract_block("kinetics", 0, blob_vector, &kinetics_vector) == FFEA_ERROR)
                return FFEA_ERROR;

            // Get map info if necessary
            if (params.num_conformations[i] > 1)
            {

                // Get map data
                if (systemreader->extract_block("maps", 0, kinetics_vector, &map_vector) == FFEA_ERROR)
                    return FFEA_ERROR;

                // Parse map data
                for (it = map_vector.begin(); it != map_vector.end(); ++it)
                {
                    systemreader->parse_map_tag(*it, map_indices, &map_fname);
                    blob_conf[i].maps.push_back(map_fname);
                    blob_conf[i].maps_conf_index_from.push_back(map_indices[0]);
                    blob_conf[i].maps_conf_index_to.push_back(map_indices[1]);
                }

                // Error check
                if (blob_conf[i].maps.size() != params.num_conformations[i] * (params.num_conformations[i] - 1))
                {
                    FFEA_error_text();
                    cout << "In blob " << i << ", expected " << params.num_conformations[i] * (params.num_conformations[i] - 1) << " maps to describe all possible switches.\n Read " << blob_conf[i].maps.size() << " maps." << endl;
                    return FFEA_ERROR;
                }
            }

            // Then, states and rates data if necessary
            if (params.num_states[i] > 1)
            {
                for (it = kinetics_vector.begin(); it != kinetics_vector.end(); ++it)
                {
                    systemreader->parse_tag(*it, lrvalue);
                    if (lrvalue[0] == "maps" || lrvalue[0] == "/maps")
                    {
                        continue;
                    }
                    else if (lrvalue[0] == "states")
                    {
                        blob_conf[i].states = lrvalue[1];
                        set_states = 1;
                    }
                    else if (lrvalue[0] == "rates")
                    {
                        blob_conf[i].rates = lrvalue[1];
                        set_rates = 1;
                    }
                }

                // Error check
                if (set_states == 0)
                {
                    FFEA_ERROR_MESSG("Expected a .states files for blob %d, as num_states = %d\n", i, params.num_states[i])
                }
                if (set_rates == 0)
                {
                    FFEA_ERROR_MESSG("Expected a .rates files for blob %d, as num_states = %d\n", i, params.num_states[i])
                }
            }
            else
            {
                blob_conf[i].states = "";
                blob_conf[i].rates = "";
            }
        }

        // Finally, get the extra blob data (solver, scale, centroid etc)
        int rotation_type = -1;

        for (it = blob_vector.begin(); it != blob_vector.end(); ++it)
        {
            systemreader->parse_tag(*it, lrvalue);

            if (lrvalue[0] == "conformation" || lrvalue[0] == "/conformation" || lrvalue[0] == "kinetics" || lrvalue[0] == "/kinetics" || lrvalue[0] == "maps" || lrvalue[0] == "/maps")
            {
                continue;
            }
            else if (lrvalue[0] == "solver")
            {
                if (lrvalue[1] == "CG")
                {
                    solver = FFEA_ITERATIVE_SOLVER;
                    mass_in_system = true;
                }
                else if (lrvalue[1] == "CG_nomass")
                {
                    solver = FFEA_NOMASS_CG_SOLVER;
                }
                else if (lrvalue[1] == "direct")
                {
                    solver = FFEA_DIRECT_SOLVER;
                    mass_in_system = true;
                }
                else if (lrvalue[1] == "masslumped")
                {
                    solver = FFEA_MASSLUMPED_SOLVER;
                    mass_in_system = true;
                }
                else
                {
                    FFEA_error_text();
                    cout << "In blob " << i << ", unrecognised solver type.\nRecognised solvers:CG, CG_nomass, direct, masslumped." << endl;
                    return FFEA_ERROR;
                }
                set_solver = 1;
            }
            else if (lrvalue[0] == "scale")
            {
                scale = atof(lrvalue[1].c_str());
                set_scale = 1;
                scale /= mesoDimensions::length;
            }
            else if (lrvalue[0] == "calc_compress")
            {
                calc_compress = atof(lrvalue[1].c_str());
            }
            else if (lrvalue[0] == "compress")
            {
                compress = atof(lrvalue[1].c_str());
            }
            else if (lrvalue[0] == "centroid" || lrvalue[0] == "centroid_pos")
            {
                blob_conf[i].set_centroid = 1;
                lrvalue[1] = boost::erase_last_copy(boost::erase_first_copy(lrvalue[1], "("), ")");
                boost::trim(lrvalue[1]);
                systemreader->split_string(lrvalue[1], blob_conf[i].centroid, ",", 3);
            }
            else if (lrvalue[0] == "velocity")
            {
                blob_conf[i].set_velocity = 1;
                lrvalue[1] = boost::erase_last_copy(boost::erase_first_copy(lrvalue[1], "("), ")");
                boost::trim(lrvalue[1]);
                systemreader->split_string(lrvalue[1], blob_conf[i].velocity, ",", 3);
                for (int ivt = 0; ivt < 3; ivt++)
                {
                    blob_conf[i].velocity[ivt] /= mesoDimensions::velocity;
                }
            }
            else if (lrvalue[0] == "rotation")
            {
                blob_conf[i].set_rotation = 1;
                lrvalue[1] = boost::erase_last_copy(boost::erase_first_copy(lrvalue[1], "("), ")");
                boost::trim(lrvalue[1]);
                if (systemreader->split_string(lrvalue[1], blob_conf[i].rotation, ",", 9) == 3)
                {
                    blob_conf[i].rotation_type = 0;
                }
                else
                {
                    blob_conf[i].rotation_type = 1;
                }
            }
        }

        // Error checking
        if (set_solver == 0)
        {
            cout << "\tBlob " << i << ", solver not set. Defaulting to CG_nomass." << endl;
            solver = FFEA_NOMASS_CG_SOLVER;
            set_solver = 1;
        }

        if (set_scale == 0)
        {
            FFEA_error_text();
            cout << "Blob " << i << ", scale not set." << endl;
            return FFEA_ERROR;
        }

        // Build blob
        // Build conformations (structural data)
        // vector3 *cent = new vector3;
        for (j = 0; j < params.num_conformations[i]; ++j)
        {
            cout << "\tConfiguring blob " << i << " conformation " << j << "..." << endl;

            if (blob_array[i][j].config(i, j, nodes[j], topology[j], surface[j],
                                        material[j], stokes[j], ssint[j], pin[j], binding[j],
                                        beads[j], scale, calc_compress, compress, solver,
                                        motion_state[j], &params, &pc_params, &ssint_matrix,
                                        &binding_matrix, rng) == FFEA_ERROR)
            {
                FFEA_ERROR_MESSG("\tError when trying to pre-initialise Blob %d, conformation %d.\n", i, j);
            }
        }

        // Clear blob vector and other vectors for next round
        motion_state.clear();
        nodes.clear();
        topology.clear();
        surface.clear();
        material.clear();
        stokes.clear();
        ssint.clear();
        binding.clear();
        pin.clear();
        beads.clear();
        scale = 1;
        solver = FFEA_NOMASS_CG_SOLVER;
        map_vector.clear();      // container for the input <map> block
        kinetics_vector.clear(); // container for the input <kinetics> block
        blob_vector.clear();     // container for the input <blob> block
        set_scale = 0;
        set_solver = 0;
        set_rates = 0;  // aux reader flag
        set_states = 0; // aux reader flag
    }

    // Blobs are now configured. Initialisation will allocate memory,
    //    and thus it may be performance wise to initialise things in the
    //    thread they will be. Hopefully will work, though that
    //    leaves us with schedule static.
    int fatal_errors = 0;
#ifdef FFEA_PARALLEL_PER_BLOB
#pragma omp parallel for default(none) schedule(static) private(i, j) reduction(+ \
                                                                                : fatal_errors) shared(blob_conf) // shared(params, blob_array, systemreader, blob_conf)

#endif
    for (i = 0; i < params.num_blobs; ++i)
    {
        for (j = 0; j < params.num_conformations[i]; ++j)
        {

            if (blob_array[i][j].init() == FFEA_ERROR)
            {
                printf("Error when initialising Blob %d, conformation %d", i, j);
                fatal_errors += 1;
            }

            // If not an active conforamtion, move to random area in infinity so ssint and stuff are not active (face linked list is not set up for deleting elements)
            if (j > 0)
            {
                blob_array[i][j].position(blob_array[i][j].get_RandU01() * 1e10, blob_array[i][j].get_RandU01() * 1e10, blob_array[i][j].get_RandU01() * 1e10);
            }
            else
            {

                // Activate all faces
                blob_array[i][j].kinetically_set_faces(true);

                // if centroid position is set, position the blob's centroid at that position. If ssint is set, move to center of box
                if (blob_conf[i].set_centroid)
                {

                    blob_conf[i].centroid[0] *= blob_array[i][j].get_scale();
                    blob_conf[i].centroid[1] *= blob_array[i][j].get_scale();
                    blob_conf[i].centroid[2] *= blob_array[i][j].get_scale();
                    vector3 dv = blob_array[i][j].position(blob_conf[i].centroid[0],
                                                           blob_conf[i].centroid[1], blob_conf[i].centroid[2]);

                    // if Blob has a number of beads, transform them too:
                    if (blob_array[i][j].get_num_beads() > 0)
                        blob_array[i][j].position_beads(dv.x, dv.y, dv.z);

                    // transform the rod, too
                    // note to future FFEA authors: yes, you have to translate your stuff here as well
                    //float shift_rod[3] = {(float)dv.x, (float)dv.y, (float)dv.z};
                    //for (i = 0; i < params.num_rods; i++) {
                    //    rod_array[i]->translate_rod(rod_array[i]->current_r, shift_rod);
                    //    rod_array[i]->translate_rod(rod_array[i]->equil_r, shift_rod);
                    //}
                }

                if (blob_conf[i].set_rotation)
                {
                    if (blob_conf[i].rotation_type == 0)
                    {
                        blob_array[i][j].rotate(blob_conf[i].rotation[0], blob_conf[i].rotation[1], blob_conf[i].rotation[2], blob_array[i][j].is_using_beads());
                    }
                    else
                    {
                        blob_array[i][j].rotate(blob_conf[i].rotation[0], blob_conf[i].rotation[1], blob_conf[i].rotation[2], blob_conf[i].rotation[3], blob_conf[i].rotation[4], blob_conf[i].rotation[5], blob_conf[i].rotation[6], blob_conf[i].rotation[7], blob_conf[i].rotation[8], blob_array[i][j].is_using_beads());
                    }
                }

                if (blob_conf[i].set_velocity)
                {
                    blob_array[i][j].velocity_all(blob_conf[i].velocity[0], blob_conf[i].velocity[1], blob_conf[i].velocity[2]);
                }

                // Set up extra nodes if necessary (STATIC structures automatically load no topology; means no internal nodes!)
                if (blob_array[i][j].get_motion_state() == FFEA_BLOB_IS_STATIC && (params.calc_steric == 1))
                {
                    blob_array[i][j].add_steric_nodes();
                }

                // set the current node positions as pos_0 for this blob, so that all rmsd values
                // are calculated relative to this conformation centred at this point in space.
                //blob_array[i][j].set_pos_0();
            }
        }

// Build kinetic system (at world level, for future potential global kinetics)
/* */
#pragma omp critical // maybe not needed, but     \
                     //   I haven't written this, \
                     //   it's still experimental, and don't have time to read.
        if (params.calc_kinetics == 1)
        {
            printf("\tInitialising kinetics for blob %d...\n", i);

            // Will load a default state if params.num_states[i] == 1
            printf("\t\tLoading kinetic states...");
            if (load_kinetic_states(blob_conf[i].states, i) == FFEA_ERROR)
            {
                FFEA_error_text();
                printf("\nProblem reading kinetic states in 'read_kinetic_states' function");
                fatal_errors += 1;
            }
            printf("...done!");

            printf("\t\tLoading kinetic rates...");
            if (load_kinetic_rates(blob_conf[i].rates, i) == FFEA_ERROR)
            {
                FFEA_error_text();
                printf("\nProblem reading kinetic rates in 'read_kinetic_rates' function");
                fatal_errors += 1;
            }
            printf("...done!\n");

            // Maps only if num_conforamtions for this blob > 1
            if (params.num_conformations[i] > 1)
            {
                printf("\t\tLoading kinetic maps...\n");
                if (load_kinetic_maps(blob_conf[i].maps, blob_conf[i].maps_conf_index_from, blob_conf[i].maps_conf_index_to, i) == FFEA_ERROR)
                {
                    FFEA_error_text();
                    printf("\nProblem reading kinetic maps in 'read_kinetic_maps' function");
                    fatal_errors += 1;
                }
                printf("\t\t...done!\n");

                printf("\t\tBuilding 'identity' maps for energy comparison...");
                if (build_kinetic_identity_maps() == FFEA_ERROR)
                {
                    FFEA_error_text();
                    printf("\nProblem reading kinetic maps in 'build_kinetic_identity_maps' function");
                    fatal_errors += 1;
                }
                printf("...done!\n");
            }
            printf("\t...done!");
        }
    }
    if (fatal_errors > 0)
    {
        cout << "There were fatal errors initialising the blobs" << endl;
        return FFEA_ERROR;
    }

    // Create rods
    rod_array = new rod::Rod *[params.num_rods];
    for (int i = 0; i < params.num_rods; i++)
    {
        // Init rod object
        rod_array[i] = rod_from_block(script_vector, i, systemreader);

        // Carry over params from script file
        rod_array[i]->viscosity = params.stokes_visc;
        rod_array[i]->timestep = params.dt;
        rod_array[i]->kT = params.kT;
        rod_array[i]->calc_noise = params.calc_noise;
        rod_array[i]->calc_steric = params.calc_steric_rod;
        rod_array[i]->calc_vdw = params.calc_vdw_rod;
        rod_array[i]->pbc = params.pbc_rod;
        rod_array[i]->flow_profile = params.flow_profile;
        for (int j=0;j<3;j++)
            rod_array[i]->flow_velocity[j] = params.flow_velocity[j];
        rod_array[i]->shear_rate = params.shear_rate;

    }

    // Create rod-blob interfaces
    rod_blob_interface_array = new rod::Rod_blob_interface *[params.num_interfaces];
    for (int i = 0; i < params.num_interfaces; i++)
    {
        rod_blob_interface_array[i] = rod_blob_interface_from_block(script_vector, i, systemreader, rod_array, blob_array);
    }

    // Set up rod-blob interfaces (in order!)
    // Iterate over 'order' variable to set up resolve order of connection positioning
    for (int order = 0; order < 999; order++)
    {
        for (int i = 0; i < params.num_interfaces; i++)
        {
            if (rod_blob_interface_array[i]->order == order)
            {
                std::cout << "Ordering interface " << i << "...\n";
                rod_blob_interface_array[i]->update_internal_state(true, true);
                if (rod_blob_interface_array[i]->ends_at_rod)
                {
                    std::cout << "Positioning rod from blob...\n";
                    rod_blob_interface_array[i]->position_rod_from_blob(false);
                    rod_blob_interface_array[i]->position_rod_from_blob(true);
                }
                else
                {
                    std::cout << "Positioning blob from rod...\n";
                    rod_blob_interface_array[i]->position_blob_from_rod();
                }
            }
        }
    }

    for (int i = 0; i < params.num_interfaces; i++)
    {
        rod_blob_interface_array[i]->set_initial_values();
    }

    // Finally, get springs
    if (params.calc_springs == 1)
    {
        if (load_springs(params.springs_fname.c_str()) != 0)
        {
            //if(load_springs(auxpath.string().c_str()) != 0) {
            //if(load_springs(lrvalue[1].c_str()) != 0) {
            FFEA_error_text();
            cout << "Problem loading springs from " << lrvalue[1] << "." << endl;
            return FFEA_ERROR;
        }
    }

    return FFEA_OK;
}

/**
 * @brief Loads the maps from a given blob
 * @param[in] vector<string> map_fnames A vector of maps to load
 * @param[in] vector<int> map_from Which conformation the maps are from
 * @param[in] vector<int> map_to Which conformation the maps go to
 * @param[in] int blob_index Which blob the maps belong to
 * @details This function reads in the maps required for kinetic switching
 * between conformations. Does no error checking for correct number of maps
 * */

int World::load_kinetic_maps(vector<string> map_fnames, vector<int> map_from, vector<int> map_to, int blob_index)
{

    unsigned int i, j, num_rows, num_cols, num_entries;
    string buf_string;
    vector<string> string_vec;
    for (i = 0; i < map_fnames.size(); ++i)
    {

        cout << "\t\t\tReading map " << i << ", " << map_fnames.at(i) << ": from conformation " << map_from.at(i) << " to " << map_to.at(i) << endl;
        ifstream fin;
        fin.open(map_fnames.at(i).c_str());
        if (fin.is_open() == false)
        {
            cout << "File " << map_fnames.at(i) << " not found. Plaese supply valid map file with correct path." << endl;
            return FFEA_ERROR;
        }

        // Check if sparse or dense
        getline(fin, buf_string);
        if (buf_string != "FFEA Kinetic Conformation Mapping File (Sparse)")
        {
            FFEA_error_text();
            cout << "In " << map_fnames.at(i) << ", expected 'FFEA Kinetic Conformation Mapping File (Sparse)'" << endl;
            cout << "but got " << buf_string << endl;
            return FFEA_ERROR;
        }

        // Get nodes to and from and check against the structures
        getline(fin, buf_string);
        boost::split(string_vec, buf_string, boost::is_space(), boost::token_compress_on);
        num_cols = atoi(string_vec.at(string_vec.size() - 1).c_str());

        if (num_cols != blob_array[blob_index][map_from.at(i)].get_num_nodes())
        {
            FFEA_error_text();
            cout << "In " << map_fnames.at(i) << ", 'num_nodes_from', " << num_cols << ", does not correspond to the number of nodes in blob " << i << " conformation " << map_from.at(i) << endl;
            return FFEA_ERROR;
        }

        getline(fin, buf_string);
        boost::split(string_vec, buf_string, boost::is_space(), boost::token_compress_on);
        num_rows = atoi(string_vec.at(string_vec.size() - 1).c_str());

        if (num_rows != blob_array[blob_index][map_to.at(i)].get_num_nodes())
        {
            FFEA_error_text();
            cout << "In " << map_fnames.at(i) << ", 'num_nodes_to', " << num_rows << ", does not correspond to the number of nodes in blob " << i << " conformation " << map_to.at(i) << endl;
            return FFEA_ERROR;
        }

        // Get num_entries
        getline(fin, buf_string);
        boost::split(string_vec, buf_string, boost::is_space(), boost::token_compress_on);
        num_entries = atoi(string_vec.at(string_vec.size() - 1).c_str());

        // Create some memory for the matrix stuff
        scalar *entries = new scalar[num_entries];
        int *key = new int[num_rows + 1];
        int *col_index = new int[num_entries];

        // Get 'map:'
        getline(fin, buf_string);

        // Read matrix
        // 'entries -'
        fin >> buf_string;
        fin >> buf_string;
        for (j = 0; j < num_entries; j++)
        {
            fin >> buf_string;
            entries[j] = atof(buf_string.c_str());
        }

        // 'key -'
        fin >> buf_string;
        fin >> buf_string;
        for (j = 0; j < num_rows + 1; j++)
        {
            fin >> buf_string;
            key[j] = atoi(buf_string.c_str());
        }

        // 'columns -'
        fin >> buf_string;
        fin >> buf_string;
        for (j = 0; j < num_entries; j++)
        {
            fin >> buf_string;
            col_index[j] = atoi(buf_string.c_str());
        }

        // Close file
        fin.close();

        // Create sparse matrix
        kinetic_map[blob_index][map_from[i]][map_to[i]].init(num_rows, num_entries, entries, key, col_index);
    }

    return FFEA_OK;
}

/**
 * @brief Builds additional maps for energy calculations
 * @details This function uses the existing kinetic maps to build
 * 'return_maps', which are used for energy calculations and comparisons
 * */

int World::build_kinetic_identity_maps()
{

    int i, j, k;

    // For each blob, build map_ij*map_ji and map_ji*map_ij so we can compare energies using only the conserved modes. Well clever this, Oliver Harlen's idea.
    // He didn't write this though!

    for (i = 0; i < params.num_blobs; ++i)
    {
        for (j = 0; j < params.num_conformations[i]; ++j)
        {
            for (k = j + 1; k < params.num_conformations[i]; ++k)
            {

                kinetic_return_map[i][j][k] = kinetic_map[i][k][j].apply(&kinetic_map[i][j][k]);
                kinetic_return_map[i][k][j] = kinetic_map[i][j][k].apply(&kinetic_map[i][k][j]);
            }
        }
    }
    return FFEA_OK;
}

/**
 * @brief Calculates kinetic rates based upon the current state of the blob
 * @details This function alters the given kinetic_rates using the energy of the system.
 * The average rate throughout the simulation should still be the given values.
 * */

int World::calculate_kinetic_rates()
{

    int i, j;
    int current_state;
    float prob_sum;

    int base_bsindex, target_bsindex, other_blob_index;
    int base_type, target_type;
    BindingSite *base_site, *target_site;

    // For each blob
    for (i = 0; i < params.num_blobs; ++i)
    {

        // Get current state
        current_state = active_blob_array[i]->get_state_index();
        //cout << "Current State = " << current_state << endl;
        // Set total probability to 0
        prob_sum = 0.0;

        // And for each state we could switch to
        for (j = 0; j < params.num_states[i]; ++j)
        {

            // No need to check if j == current_state
            if (current_state == j)
            {
                continue;
            }

            // Or if the base rate is zero
            if (kinetic_base_rate[i][current_state][j] == 0)
            {
                kinetic_rate[i][current_state][j] = kinetic_base_rate[i][current_state][j];
                continue;
            }

            // What type of state change do we have? (these options should be exclusive. make sure of this in initialisation)
            if (kinetic_state[i][current_state].get_conformation_index() != kinetic_state[i][j].get_conformation_index())
            {

                // Conformation change! Kinetic switch is dependent upon the energies (or they will be at least!)
                kinetic_rate[i][current_state][j] = kinetic_base_rate[i][current_state][j];
            }
            else if (!kinetic_state[i][current_state].is_bound() && kinetic_state[i][j].is_bound())
            {

                // Binding event! Kinetic switch is constant but a step function dependent upon distance from the potential binding sites. Entropy taken into account by simulation
                // Initialise to zero in case of no sites in range
                kinetic_rate[i][current_state][j] = 0.0;

                // Get the base and target types
                base_type = kinetic_state[i][j].get_base_bsite_type();
                target_type = kinetic_state[i][j].get_target_bsite_type();

                // Scan all sites on this blob
                for (base_bsindex = 0; base_bsindex < active_blob_array[i]->num_binding_sites; ++base_bsindex)
                {
                    base_site = active_blob_array[i]->get_binding_site(base_bsindex);

                    // If wrong type, move on
                    if (base_site->get_type() != base_type)
                    {
                        continue;
                    }

                    // Else, scan all other binding sites in the world
                    for (other_blob_index = 0; other_blob_index < params.num_blobs; ++other_blob_index)
                    {

                        // If same blob, continue (for now)
                        if (i == other_blob_index)
                        {
                            continue;
                        }

                        // Scan all sites on this blob too
                        for (target_bsindex = 0; target_bsindex < active_blob_array[other_blob_index]->num_binding_sites; ++target_bsindex)
                        {

                            target_site = active_blob_array[other_blob_index]->get_binding_site(target_bsindex);

                            // If wrong type, move on
                            if (target_site->get_type() != target_type)
                            {
                                continue;
                            }

                            // We've got 2 compatible sites! Are they in range?
                            if (BindingSite::sites_in_range(*base_site, *target_site))
                            {

                                // Success! Set rates and bsites into the states
                                kinetic_rate[i][current_state][j] = kinetic_base_rate[i][current_state][j];
                                kinetic_state[i][j].set_sites(base_site, target_site);

                                // And return from this crazy loop
                                other_blob_index = params.num_blobs;
                                base_bsindex = active_blob_array[i]->num_binding_sites;
                                break;
                            }
                        }
                    }
                }
            }
            else if (kinetic_state[i][current_state].is_bound() && !kinetic_state[i][j].is_bound())
            {

                // Unbnding event! Kinetic switch is constant
                kinetic_rate[i][current_state][j] = kinetic_base_rate[i][current_state][j];
            }
            else
            {

                // Identity event. Nothing changes here either
                kinetic_rate[i][current_state][j] = kinetic_base_rate[i][current_state][j];
            }

            prob_sum += kinetic_rate[i][current_state][j];
        }

        // Finally, the probability of staying put
        if (prob_sum > 1.0)
        {
            FFEA_ERROR_MESSG("Although your original switching probabilities for blob %d totalled < 1.0, after rescaling they have gone > 1.0. Lower your 'kinetic_update' parameter!", i)
        }

        kinetic_rate[i][current_state][current_state] = 1 - prob_sum;
    }

    return FFEA_OK;
}

/**
 * @brief Selects a new states based on the current kinetic rates
 * @param[in] int blob_index Which blob wants to switch
 * @param[in] int *target A list of potential states
 * @details This function randomly chooses a state to switch to
 * out of the given allowed taget states based upon the current kinetic
 * rates.
 * */

int World::choose_new_kinetic_state(int blob_index, int *target)
{

    int i;
    scalar switch_check;

    // Make some bins
    vector<scalar[2]> bin(params.num_states[blob_index]);
    scalar total = 0.0;
    for (i = 0; i < params.num_states[blob_index]; ++i)
    {

        // Lower limit
        bin[i][0] = total;

        // Upper limit
        total += kinetic_rate[blob_index][active_blob_array[blob_index]->get_state_index()][i];
        bin[i][1] = total;
    }

    // Round up in case of numerical problems
    bin[params.num_states[blob_index] - 1][1] = 1.0;

    // Get a random number
    switch_check = kinetic_rng->RandU01();

    // See which bin this falls into
    for (i = 0; i < params.num_states[blob_index]; ++i)
    {
        if (switch_check > bin[i][0] && switch_check <= bin[i][1])
        {
            *target = i;
        }
    }
    return FFEA_OK;
}

int World::load_kinetic_states(string states_fname, int blob_index)
{

    int i, j, num_states, conf_index, from, to; //, site_index;
    const int MAX_BUF_SIZE = 255;
    char buf[MAX_BUF_SIZE];
    string buf_string;
    vector<string> sline;
    vector<string>::iterator it;

    // Load a default, single state
    if (states_fname == "")
    {

        kinetic_state[blob_index] = new KineticState[1];
        kinetic_state[blob_index][0].init();

        return FFEA_OK;
    }

    // Else

    // Open the file
    ifstream fin;
    fin.open(states_fname);
    if (fin.fail())
    {
        FFEA_ERROR_MESSG("'states_fname' %s not found\n", states_fname.c_str())
    }

    cout << "\n\t\tReading in Kinetic States file: " << states_fname << endl;

    // Get header stuff and check for errors
    fin.getline(buf, MAX_BUF_SIZE);

    if (strcmp(buf, "ffea kinetic states file") != 0)
    {
        FFEA_ERROR_MESSG("\nExpected 'ffea kinetic states file' as first line. This may not be an FFEA kinetic states file\n")
    }

    // Get num_states
    getline(fin, buf_string);
    boost::trim(buf_string);
    boost::split(sline, buf_string, boost::is_space(), boost::token_compress_on);
    ;
    num_states = atoi(sline.at(1).c_str());

    if (num_states != params.num_states[blob_index])
    {
        FFEA_ERROR_MESSG("\nnum_states defined in '%s', %d, does not correspond to the initial script file, %d.\n", states_fname.c_str(), num_states, params.num_states[blob_index])
    }

    // Get 'states:'
    getline(fin, buf_string);

    // Create state objects
    kinetic_state[blob_index] = new KineticState[num_states];

    // Get actual states (each line varies)
    for (i = 0; i < num_states; ++i)
    {

        // Get conformation_index first
        getline(fin, buf_string);
        boost::trim(buf_string);
        boost::split(sline, buf_string, boost::is_space(), boost::token_compress_on);

        conf_index = atoi(sline.at(1).c_str());

        if (conf_index < 0 || conf_index >= params.num_conformations[blob_index])
        {
            FFEA_ERROR_MESSG("In %s, state %d, conf_index is out of range (0 < conf_index < %d).\n", states_fname.c_str(), i, params.num_conformations[blob_index])
        }

        // Now get bound binding sites
        sline.clear();
        getline(fin, buf_string);
        boost::trim(buf_string);
        boost::split(sline, buf_string, boost::is_space(), boost::token_compress_on);

        // Check line consistency

        // No sites defined
        if (sline.size() == 1)
        {
            from = -1;
            to = -1;
        }
        else
        {

            if (sline.size() - 1 != 2)
            {
                FFEA_ERROR_MESSG("In %s, state %d, binding line format error. Should be 'binding from_index to_index\n", states_fname.c_str(), i)
            }

            sline.erase(sline.begin());

            // Reset counter
            j = -1;
            for (it = sline.begin(); it != sline.end(); ++it)
            {
                j++;
                if (j == 0)
                {
                    from = atoi((*it).c_str());
                }
                else if (j == 1)
                {
                    to = atoi((*it).c_str());

                    if (from >= binding_matrix.num_interaction_types || to >= binding_matrix.num_interaction_types)
                    {
                        FFEA_ERROR_MESSG("In %s, state %d, binding from type %d or to type %d is > num_interaction_types %d (check binding_matrix).\n", states_fname.c_str(), i, from, to, binding_matrix.num_interaction_types)
                    }

                    if (binding_matrix.interaction[from][to] != true)
                    {
                        FFEA_ERROR_MESSG("In %s, state %d, binding from type %d to type %d is not allowed (check binding_matrix).\n", states_fname.c_str(), i, from, to)
                    }
                }
            }
        }

        // Initialise kinetic state from this
        kinetic_state[blob_index][i].init(conf_index, from, to);
    }

    // Close and return
    fin.close();
    return FFEA_OK;
}

int World::load_kinetic_rates(string rates_fname, int blob_index)
{

    int i, j, num_states;
    char *crap;
    char buf[255];
    string buf_string;
    vector<string> sline;
    vector<string>::iterator it;
    FILE *fin;

    // Load a default, single rate
    if (rates_fname == "")
    {

        // Create rates matrix
        kinetic_rate[blob_index] = new scalar *[1];
        kinetic_base_rate[blob_index] = new scalar *[1];

        kinetic_rate[blob_index][0] = new scalar[1];
        kinetic_base_rate[blob_index][0] = new scalar[1];

        kinetic_rate[blob_index][0][0] = 1.0;
        kinetic_base_rate[blob_index][0][0] = 1.0;

        return FFEA_OK;
    }

    // Open the file
    fin = fopen(rates_fname.c_str(), "rb");

    // Get header stuff and check for errors
    crap = fgets(buf, 255, fin);
    if (strcmp(buf, "ffea kinetic rates file\n") != 0)
    {
        FFEA_ERROR_MESSG("\nExpected 'ffea kinetic rates file' as first line. This may not be an FFEA kinetic rates file\n")
    }

    if (fscanf(fin, "num_states %d\n", &num_states) != 1)
    {
        FFEA_ERROR_MESSG("\nExpected 'num_states %%d' as second line. Unable to read further.\n")
    }
    if (num_states != params.num_states[blob_index])
    {
        FFEA_ERROR_MESSG("\nnum_states defined in '%s', %d, does not correspond to the initial script file, %d.\n", rates_fname.c_str(), num_states, params.num_states[blob_index])
    }
    crap = fgets(buf, 255, fin);

    // Create rates matrix
    kinetic_rate[blob_index] = new scalar *[num_states];
    kinetic_base_rate[blob_index] = new scalar *[num_states];
    for (i = 0; i < num_states; ++i)
    {
        kinetic_rate[blob_index][i] = new scalar[num_states];
        kinetic_base_rate[blob_index][i] = new scalar[num_states];
    }

    scalar total_prob;

    // Get each state's rates and check total probability is conserved
    for (i = 0; i < num_states; ++i)
    {
        total_prob = 0.0;

        // Get a line and split it
        crap = fgets(buf, 255, fin);
        boost::split(sline, buf, boost::is_any_of(" "), boost::token_compress_on);
        if (sline.size() > num_states)
        {
            FFEA_ERROR_MESSG("\nState %d contains %zd rate values, instead of 'num_states', %d.\n", i, sline.size(), num_states)
        }

        j = -1;
        for (it = sline.begin(); it != sline.end(); ++it)
        {

            // Increment counter
            j++;
            kinetic_base_rate[blob_index][i][j] = atof((*it).c_str());
            kinetic_base_rate[blob_index][i][j] *= mesoDimensions::time;

            // Change to probabilities and ignore diagonal
            kinetic_base_rate[blob_index][i][j] *= params.dt * params.kinetics_update;
            if (i != j)
            {
                total_prob += kinetic_base_rate[blob_index][i][j];
            }
        }

        // Prob of not switching (for completion)
        if (total_prob > 1)
        {
            FFEA_error_text();
            cout << "P(switch_state in kinetic update period) = rate(switch_state)(Hz) * dt * kinetics_update" << endl;
            cout << "Due to the size of your rates, your timestep, and your kinetic_update value, the total probability of changing states each kinetic update period is greater than one." << endl;
            cout << "Best solution - Reduce 'kinetics_update' parameter" << endl;
            return FFEA_ERROR;
        }
        kinetic_base_rate[blob_index][i][i] = 1 - total_prob;
    }
    return FFEA_OK;
}

void World::print_kinetic_rates_to_screen(int type)
{

    int i, j, k;
    cout << "Kinetic Rates:" << endl;
    for (i = 0; i < params.num_blobs; ++i)
    {
        cout << "\tBlob " << i << ":\n"
             << endl;
        cout << "\tto";
        for (j = 0; j < params.num_states[i]; ++j)
        {
            cout << "\t" << j;
        }
        cout << endl
             << "from" << endl;
        for (j = 0; j < params.num_states[i]; ++j)
        {
            cout << j << "\t\t";
            for (k = 0; k < params.num_states[i]; ++k)
            {
                if (type == 0)
                {
                    cout << kinetic_base_rate[i][j][k] << "\t";
                }
                else if (type == 1)
                {
                    cout << kinetic_rate[i][j][k] << "\t";
                }
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
}

/* */
void World::get_system_CoM(vector3 *system_CoM)
{
    system_CoM->x = 0;
    system_CoM->y = 0;
    system_CoM->z = 0;
    scalar total_mass = 0;
    for (int i = 0; i < params.num_blobs; i++)
    {
        vector3 com;
        active_blob_array[i]->get_CoM(&com);
        system_CoM->x += com.x * active_blob_array[i]->get_mass();
        system_CoM->y += com.y * active_blob_array[i]->get_mass();
        system_CoM->z += com.z * active_blob_array[i]->get_mass();

        total_mass += active_blob_array[i]->get_mass();
    }
    system_CoM->x /= total_mass;
    system_CoM->y /= total_mass;
    system_CoM->z /= total_mass;
}

/* */
void World::get_system_centroid(vector3 *centroid)
{
    /** Blob centroid */
    centroid->x = 0;
    centroid->y = 0;
    centroid->z = 0;
    int total_num_nodes = 0;
    for (int i = 0; i < params.num_blobs; i++)
    {
        vector3 cen;
        active_blob_array[i]->get_centroid(&cen);
        centroid->x += cen.x * active_blob_array[i]->get_num_nodes();
        centroid->y += cen.y * active_blob_array[i]->get_num_nodes();
        centroid->z += cen.z * active_blob_array[i]->get_num_nodes();

        total_num_nodes += active_blob_array[i]->get_num_nodes();
    }

    /** Rod centroid */
    for (int i = 0; i < params.num_rods; i++)
    {
        float rod_centroid[3];
        rod_array[i]->get_centroid(rod_array[i]->current_r, rod_centroid);
        /** I'm leaving it like this and there's nothing you can do about it */
        centroid->x += rod_centroid[0] * rod_array[i]->num_nodes;
        centroid->y += rod_centroid[1] * rod_array[i]->num_nodes;
        centroid->z += rod_centroid[2] * rod_array[i]->num_nodes;
        total_num_nodes += rod_array[i]->num_nodes;
    }

    centroid->x /= total_num_nodes;
    centroid->y /= total_num_nodes;
    centroid->z /= total_num_nodes;
}

void World::get_system_dimensions(vector3 *dimension)
{
    dimension->x = 0;
    dimension->y = 0;
    dimension->z = 0;

    vector3 min, max;
    min.x = INFINITY;
    min.y = INFINITY;
    min.z = INFINITY;
    max.x = -1 * INFINITY;
    max.y = -1 * INFINITY;
    max.z = -1 * INFINITY;

    vector3 blob_min, blob_max;
    for (int i = 0; i < params.num_blobs; i++)
    {
        active_blob_array[i]->get_min_max(&blob_min, &blob_max);
        if (blob_min.x < min.x)
        {
            min.x = blob_min.x;
        }
        if (blob_min.y < min.y)
        {
            min.y = blob_min.y;
        }
        if (blob_min.z < min.z)
        {
            min.z = blob_min.z;
        }

        if (blob_max.x > max.x)
        {
            max.x = blob_max.x;
        }

        if (blob_max.y > max.y)
        {
            max.y = blob_max.y;
        }

        if (blob_max.z > max.z)
        {
            max.z = blob_max.z;
        }
    }

    float rod_min[3], rod_max[3];
    for (int i = 0; i < params.num_rods; i++)
    {
        rod_array[i]->get_min_max(rod_array[i]->current_r, rod_min, rod_max);
        if (rod_max[0] > max.x)
        {
            max.x = rod_max[0];
        }
        if (rod_max[1] > max.y)
        {
            max.y = rod_max[1];
        }
        if (rod_max[2] > max.z)
        {
            max.z = rod_max[2];
        }
        if (rod_min[0] < min.x)
        {
            min.x = rod_min[0];
        }
        if (rod_min[1] < min.y)
        {
            min.y = rod_min[1];
        }
        if (rod_min[2] < min.z)
        {
            min.z = rod_min[2];
        }
    }

    dimension->x = max.x - min.x;
    dimension->y = max.y - min.y;
    dimension->z = max.z - min.z;
}

int World::get_num_blobs()
{

    return params.num_blobs;
}

int World::load_springs(const char *fname)
{

    int i, crap;
    FILE *in = NULL;
    const int max_line_size = 50;
    char line[max_line_size];

    // open the spring file
    if ((in = fopen(fname, "rb")) == NULL)
    {
        FFEA_FILE_ERROR_MESSG(fname)
    }
    printf("\tReading in springs file: %s\n", fname);

    // first line should be the file type "ffea springs file"
    if (fgets(line, max_line_size, in) == NULL)
    {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading first line of spring file\n")
    }
    if (strcmp(line, "ffea springs file\n") != 0)
    {
        fclose(in);
        FFEA_ERROR_MESSG("This is not a 'ffea spring file' (read '%s') \n", line)
    }

    // read in the number of springs in the file
    if (fscanf(in, "num_springs %d\n", &num_springs) != 1)
    {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading number of springs\n")
    }
    printf("\t\tNumber of springs = %d\n", num_springs);

    // Allocate memory for springs
    spring_array = new Spring[num_springs];

    // Read in next line
    crap = fscanf(in, "springs:\n");
    for (i = 0; i < num_springs; ++i)
    {
        if (fscanf(in, "%lf %lf %d %d %d %d %d %d\n", &spring_array[i].k, &spring_array[i].l, &spring_array[i].blob_index[0], &spring_array[i].blob_index[1], &spring_array[i].conformation_index[0], &spring_array[i].conformation_index[1], &spring_array[i].node_index[0], &spring_array[i].node_index[1]) != 8)
        {
            FFEA_error_text();
            printf("Problem reading spring data from %s. Format is:\n\n", fname);
            printf("ffea spring file\nnum_springs ?\n");
            printf("k l blob_index_0 blob_index 1 conformation_index_0 conformation_index_1 node_index_0 node_index_1\n\n");
            return FFEA_ERROR;
        }
        spring_array[i].k *= mesoDimensions::area / mesoDimensions::Energy;
        spring_array[i].l /= mesoDimensions::length;

        // Flag on blob to say springs are present
        active_blob_array[spring_array[i].blob_index[0]]->set_springs_on_blob(true);
        active_blob_array[spring_array[i].blob_index[1]]->set_springs_on_blob(true);

        // Error checking
        for (int j = 0; j < 2; ++j)
        {
            if (spring_array[i].blob_index[j] >= params.num_blobs || spring_array[i].blob_index[j] < 0)
            {
                FFEA_error_text();
                printf("In spring %d, blob index %d is out of bounds given the number of blobs defined (%d). Please fix. Remember, indexing starts at ZERO!\n", i, j, params.num_blobs);
                return FFEA_ERROR;
            }
            if (spring_array[i].conformation_index[j] >= params.num_conformations[spring_array[i].blob_index[j]] || spring_array[i].conformation_index[j] < 0)
            {
                FFEA_error_text();
                printf("In spring %d, conformation index %d is out of bounds given the number of conformations defined in blob %d (%d). Please fix. Remember, indexing starts at ZERO!\n", i, j, spring_array[i].blob_index[j], params.num_conformations[spring_array[i].blob_index[j]]);
                return FFEA_ERROR;
            }
            if (spring_array[i].node_index[j] >= blob_array[spring_array[i].blob_index[j]][spring_array[i].conformation_index[j]].get_num_nodes() || spring_array[i].node_index[j] < 0)
            {
                FFEA_error_text();
                printf("In spring %d, node index %d is out of bounds given the number of nodes defined in blob %d, conformation %d (%d). Please fix. Remember, indexing starts at ZERO!\n", i, j, spring_array[i].blob_index[j], spring_array[i].conformation_index[j], blob_array[spring_array[i].blob_index[j]][spring_array[i].conformation_index[j]].get_num_nodes());
                return FFEA_ERROR;
            }
            if (spring_array[i].k < 0)
            {
                FFEA_error_text();
                printf("In spring %d, spring constant, %e, < 0. This is not going to end well for you...\n", i, spring_array[i].k);
                return FFEA_ERROR;
            }
            if (spring_array[i].l < 0)
            {
                FFEA_error_text();
                printf("In spring %d, spring equilibrium length, %e, < 0. Reverse node definitions for consistency.\n", i, spring_array[i].l);
                return FFEA_ERROR;
            }
        }
        if (spring_array[i].blob_index[0] == spring_array[i].blob_index[1] && spring_array[i].conformation_index[0] == spring_array[i].conformation_index[1] && spring_array[i].node_index[0] == spring_array[i].node_index[1])
        {
            FFEA_error_text();
            printf("In spring %d, spring is connected to same node on same blob on same conformation. Will probably cause an force calculation error.\n", i);
            return FFEA_ERROR;
        }
        if (spring_array[i].blob_index[0] == spring_array[i].blob_index[1] && spring_array[i].conformation_index[0] != spring_array[i].conformation_index[1])
        {
            FFEA_error_text();
            printf("In spring %d, spring is connected two conformations of the same blob (blob_index = %d). This cannot happen as conformations are mutually exclusive.\n", i, spring_array[i].blob_index[0]);
            return FFEA_ERROR;
        }
    }

    fclose(in);

    // Inititalise the energy array (move to a solver in the future, like the ssint)
    springfieldenergy = new scalar *[params.num_blobs];
    for (i = 0; i < params.num_blobs; ++i)
    {
        springfieldenergy[i] = new scalar[params.num_blobs];
    }
    printf("\t\tRead %d springs from %s\n", num_springs, fname);
    activate_springs();
    return 0;
}

rod::Rod_blob_interface *World::rod_blob_interface_from_block(vector<string> block, int interface_id, FFEA_input_reader *systemreader, rod::Rod **rod_array, Blob **blob_array)
{

    string tag_out[2];
    bool coupling_parent = false;

    int block_no = -1;
    int coupling_counter = 0;

    bool ends_at_rod = true;
    int blob_id;
    int rod_id;
    int rod_node_id;
    int blob_element_id;
    int blob_node_ids[3];
    int from_index;
    int to_index;
    int order;
    float rotation[3];
    float node_weighting[3];

    for (auto &tag_str : block)
    {

        systemreader->parse_tag(tag_str, tag_out);

        // Are we in a <coupling> block?
        //if (tag_out[0] == "blob"){rod_parent = false;}
        if (tag_out[0] == "coupling type")
        {
            coupling_parent = true;
            block_no += 1;
        }

        if (block_no != interface_id)
        {
            continue;
        }

        // Set filename
        //if (tag_out[0] == "output" && coupling_parent){ current_rod->change_filename(tag_out[1]); }

        // parse coupling block
        if (tag_out[0] == "coupling type" && coupling_parent && (tag_out[1] == "rod-to-blob" || tag_out[1] == "blob-to-rod"))
        {
            std::cout << "Coupling data parsing...";
            if (tag_out[1] == "rod-to-blob")
            {
                ends_at_rod = false;
            }
            vector<string> sub_block;
            systemreader->extract_block("coupling", interface_id, block, &sub_block, true);
            string sub_tag_out[2];
            for (auto &sub_tag_str : sub_block)
            {
                systemreader->parse_tag(sub_tag_str, sub_tag_out);
                if ((sub_tag_out[0] == "rod_id"))
                {
                    rod_id = stof(sub_tag_out[1]);
                }
                if ((sub_tag_out[0] == "blob_id"))
                {
                    blob_id = stof(sub_tag_out[1]);
                }
                if ((sub_tag_out[0] == "rod_node_id"))
                {
                    rod_node_id = stof(sub_tag_out[1]);
                }
                if ((sub_tag_out[0] == "blob_element_id"))
                {
                    blob_element_id = stof(sub_tag_out[1]);
                }

                if (sub_tag_out[0] == "blob_node_ids")
                {
                    sub_tag_out[1] = boost::erase_last_copy(boost::erase_first_copy(sub_tag_out[1], "("), ")");
                    systemreader->split_string(sub_tag_out[1], blob_node_ids, ",", 3);
                }

                if ((sub_tag_out[0] == "order"))
                {
                    order = stoi(sub_tag_out[1]);
                }

                if (sub_tag_out[0] == "rotation")
                {
                    scalar rotation_scalar[3];
                    sub_tag_out[1] = boost::erase_last_copy(boost::erase_first_copy(sub_tag_out[1], "("), ")");
                    systemreader->split_string(sub_tag_out[1], rotation_scalar, ",", 3);
                    std::copy(rotation_scalar, rotation_scalar + 3, rotation);
                }

                if (sub_tag_out[0] == "node_weighting")
                {
                    scalar node_weighting_scalar[3];
                    sub_tag_out[1] = boost::erase_last_copy(boost::erase_first_copy(sub_tag_out[1], "("), ")");
                    systemreader->split_string(sub_tag_out[1], node_weighting_scalar, ",", 3);
                    std::cout << "node weighting = [" << node_weighting_scalar[0] << ", " << node_weighting_scalar[1] << ", " << node_weighting_scalar[2] << "\n";
                    std::copy(node_weighting_scalar, node_weighting_scalar + 3, node_weighting);
                }
            }
            coupling_counter += 1;
            break;
        }
    }

    if (ends_at_rod == false)
    {
        to_index = rod_node_id;
        from_index = blob_element_id;
    }
    else
    {
        to_index = blob_element_id;
        from_index = rod_node_id;
    }

    rod::Rod *connected_rod_ptr = rod_array[rod_id];
    Blob *connected_blob_ptr = blob_array[blob_id];

    rod::Rod_blob_interface *curr_rbi = new rod::Rod_blob_interface(connected_rod_ptr, connected_blob_ptr, ends_at_rod, to_index, from_index, blob_node_ids, rotation, node_weighting, order);

    return curr_rbi;
}

// Returns a pointer to a rod object from an .ffea script that's already
// been converted into a vector<string>. Also needs the block_id. Before
// this can be used, we need a way to get the number of rods specified
// in the file, so we can allocate an array for their pointers and know
// what block_ids to assign.
rod::Rod *World::rod_from_block(vector<string> block, int block_id, FFEA_input_reader *systemreader)
{

    // Find trajectory file
    string tag_out[2];
    string filename;
    string out_filename;
    int rod_block_no = -1; // start indexing rod blocks from 0 (we add 1 to this if rod is found)
    bool restart = false;
    for (auto &tag_str : block)
    {
        systemreader->parse_tag(tag_str, tag_out);
        if (tag_out[0] == "input")
        {
            rod_block_no += 1;
            if (rod_block_no == block_id)
            {
                filename = tag_out[1];
            }
        }
    }

    // If an output file already exists, we're on a restart, so load that instead
    rod_block_no = -1;
    for (auto &tag_str : block)
    {
        systemreader->parse_tag(tag_str, tag_out);
        if (tag_out[0] == "output")
        {
            rod_block_no += 1;
            if (rod_block_no == block_id)
            {
                out_filename = tag_out[1];
            }
        }
    }

    ifstream out_test(out_filename);
    if (out_test.good())
    {
        filename = out_filename;
        restart = true;
    }

    // Create rod object
    //rod::Rod *current_rod;
    rod::Rod *current_rod = new rod::Rod(filename, block_id);
    current_rod->load_header(filename);
    current_rod->load_contents(filename);
    current_rod->set_units();

    if (restart)
    {
        current_rod->restarting = true;
    }

    bool rod_parent = false;
    rod_block_no = -1; // start indexing rod blocks from 0

    int coupling_counter = 0;
    for (auto &tag_str : block)
    {
        systemreader->parse_tag(tag_str, tag_out);

        // Are we in a <rod> block?
        if (tag_out[0] == "blob")
        {
            rod_parent = false;
        }
        if (tag_out[0] == "rod")
        {
            rod_parent = true;
            rod_block_no += 1;
        }

        if (rod_block_no != current_rod->rod_no)
        {
            continue;
        }

        // Set filename
        if (tag_out[0] == "output" && rod_parent && !restart)
        {
            current_rod->change_filename(tag_out[1]);
        }

        // Scale rod
        if (tag_out[0] == "scale" && rod_parent && !restart)
        {
            float scale = stof(tag_out[1]);
            current_rod->scale_rod(scale);
        }

        // Move centroid
        if (tag_out[0] == "centroid_pos" && rod_parent && !restart)
        {
            // get centroid and convert it to array
            scalar centroid_pos[3];
            float converted_centroid[3];
            tag_out[1] = boost::erase_last_copy(boost::erase_first_copy(tag_out[1], "("), ")");
            systemreader->split_string(tag_out[1], centroid_pos, ",", 3);
            // convert to floats
            std::copy(centroid_pos, centroid_pos + 3, converted_centroid);
            // set centroid
            current_rod->translate_rod(current_rod->current_r, converted_centroid);
            current_rod->translate_rod(current_rod->equil_r, converted_centroid);
        }

        // Rotate rod
        if (tag_out[0] == "rotation" && rod_parent && !restart)
        {
            // get centroid and convert it to array
            scalar rotation[3];
            tag_out[1] = boost::erase_last_copy(boost::erase_first_copy(tag_out[1], "("), ")");
            systemreader->split_string(tag_out[1], rotation, ",", 3);
            // convert to floats
            float converted_rotation[3];
            std::copy(rotation, rotation + 3, converted_rotation);
            // rotate that bad boy
            current_rod->rotate_rod(converted_rotation);
        }

        // van der Waals interaction sites (.rodvdw file)
        if (tag_out[0] == "vdw" && rod_parent && !restart && params.calc_vdw_rod)
        {
            current_rod->load_vdw(tag_out[1]);
        }

    }

    return current_rod;
}

/** Populate neighbour lists of two rods, a and b.
 *
 * Loops over every element of both rods; O(N^2).
*/
void World::update_rod_steric_nbr_lists(rod::Rod *rod_a, rod::Rod *rod_b)
{
    float r_a[3] = {0};
    float r_b[3] = {0};
    float p_a[3] = {0};
    float p_b[3] = {0};

    rod_a->check_nbr_list_dim(rod_a->steric_nbrs);
    rod_b->check_nbr_list_dim(rod_b->steric_nbrs);

    if (rod::dbg_print)
        std::cout << "Updating neighbour lists of rods " << rod_a->rod_no << " and " << rod_b->rod_no << std::endl;

    for (int elem_a = 0; elem_a < rod_a->get_num_nodes() - 1; elem_a++)
    {
        for (int elem_b = 0; elem_b < rod_b->get_num_nodes() - 1; elem_b++)
        {
            rod_a->get_r(elem_a, r_a, false);
            rod_b->get_r(elem_b, r_b, false);
            rod_a->get_p(elem_a, p_a, false);
            rod_b->get_p(elem_b, p_b, false);

            // assign to both elements
            rod::set_steric_nbrs(
                rod_a->rod_no,
                rod_b->rod_no,
                elem_a,
                elem_b,
                p_a,
                p_b,
                r_a,
                r_b,
                rod_a->get_radius(elem_a),
                rod_b->get_radius(elem_b),
                rod_a->steric_nbrs.at(elem_a),
                rod_b->steric_nbrs.at(elem_b),
                params.pbc_rod,
                {(float)box_dim.x, (float)box_dim.y, (float)box_dim.z});
        }
    }
}

// ! - currently broken
void World::update_rod_vdw_nbr_lists(rod::Rod *rod_a, rod::Rod *rod_b, SSINT_matrix *lj_matrix)
{
    float r_a[3] = {0};
    float r_b[3] = {0};
    float p_a[3] = {0};
    float p_b[3] = {0};

    if (rod::dbg_print)
        std::cout << "Updating vdw neighbour lists of rods " << rod_a->rod_no << " and " << rod_b->rod_no << std::endl;

    rod_a->check_nbr_list_dim(rod_a->vdw_nbrs);
    rod_b->check_nbr_list_dim(rod_b->vdw_nbrs);

    for (auto &site_a : rod_a->vdw_sites)
    {

        for (auto &site_b : rod_b->vdw_sites)
        {
            int elem_a = site_a.elem_id;
            int elem_b = site_b.elem_id;

            site_a.update_position(rod_a->current_r);
            site_b.update_position(rod_b->current_r);

            rod_a->get_r(elem_a, r_a, false);
            rod_b->get_r(elem_b, r_b, false);
            rod_a->get_p(elem_a, p_a, false);
            rod_b->get_p(elem_b, p_b, false);

            map<string, scalar> pmap = lj_matrix->get_SSINT_params(site_a.vdw_type, site_b.vdw_type);

            rod::set_vdw_nbrs(
                site_a,
                site_b,
                p_a,
                p_b,
                r_a,
                r_b,
                rod_a->get_radius(elem_a),
                rod_b->get_radius(elem_b),
                rod_a->vdw_nbrs.at(elem_a),
                rod_b->vdw_nbrs.at(elem_b),
                params.pbc_rod,
                {(float)box_dim.x, (float)box_dim.y, (float)box_dim.z},
                params.ssint_cutoff,
                pmap["Emin"],
                pmap["Rmin"]);
        }
    }
}

// If rod centroid leaves simulation box, apply PBC wrap
void World::rod_pbc_wrap(rod::Rod *current_rod, std::vector<float> dim)
{
    for (int i = 0; i < params.num_rods; i++)
    {
        float rod_centroid[3] = {0};
        float rod_pbc_shift[3] = {0};
        current_rod->get_centroid(current_rod->current_r, rod_centroid);

        if (rod_centroid[0] < 0)
            rod_pbc_shift[0] += dim[0];
        else if (rod_centroid[0] > dim[0])
            rod_pbc_shift[0] -= dim[0];

        if (rod_centroid[1] < 0)
            rod_pbc_shift[1] += dim[1];
        else if (rod_centroid[1] > dim[1])
            rod_pbc_shift[1] -= dim[1];

        if (rod_centroid[2] < 0)
            rod_pbc_shift[2] += dim[2];
        else if (rod_centroid[2] > dim[2])
            rod_pbc_shift[2] -= dim[2];

        current_rod->translate_rod(current_rod->current_r, rod_pbc_shift);
    }
}

void World::rod_box_length_check(rod::Rod *current_rod, std::vector<float> dim)
{
    float e2e = current_rod->end_to_end_length();
    float cl = current_rod->contour_length();

    if ((e2e > 0.5 * dim[0]) || (e2e > 0.5 * dim[1]) || (e2e > 0.5 * dim[2]))
    {
        std::string msg = "box_dim: (" +
            std::to_string(dim[0]) + ", " +
            std::to_string(dim[1]) + ", " +
            std::to_string(dim[2]) + ")"
            "Rod " + std::to_string(current_rod->rod_no) +
            " end-to-end length " + std::to_string(e2e) + " exceeds half a box dimension.";
        throw std::runtime_error(msg);
    }

    if ((cl > 0.5 * dim[0]) || (cl > 0.5 * dim[1]) || (cl > 0.5 * dim[2]))
    {
        printf("\tWARNING: Rod %d contour length of %.6e exceeds half a box dimension. "
            "Interactions over boundaries may behave weirdly if it fully extends.", current_rod->rod_no, cl);
    }
}

void World::activate_springs()
{
    for (int i = 0; i < num_springs; ++i)
    {

        // If both ends of spring are active molecules, activate! This could probably be done more quickly with pointers if necessary in future
        if (spring_array[i].conformation_index[0] == active_blob_array[spring_array[i].blob_index[0]]->conformation_index && spring_array[i].conformation_index[1] == active_blob_array[spring_array[i].blob_index[1]]->conformation_index)
        {
            spring_array[i].am_i_active = true;
        }
        else
        {
            spring_array[i].am_i_active = false;
        }
    }
}

int World::apply_springs()
{
    scalar force_mag;
    vector3 n1, n0, force0, force1, sep, sep_norm;
    for (int i = 0; i < num_springs; ++i)
    {
        if (spring_array[i].am_i_active == true)
        {
            active_blob_array[spring_array[i].blob_index[1]]->get_node(spring_array[i].node_index[1], n1.data);
            active_blob_array[spring_array[i].blob_index[0]]->get_node(spring_array[i].node_index[0], n0.data);
            sep.x = n1.x - n0.x;
            sep.y = n1.y - n0.y;
            sep.z = n1.z - n0.z;

            try
            {
                arr3Normalise2<scalar, arr3>(sep.data, sep_norm.data);
            }
            catch (int e)
            {

                // If zero magnitude, we're ok
                if (e == -1)
                {
                    sep_norm.x = 0.0;
                    sep_norm.y = 0.0;
                    sep_norm.z = 0.0;
                }
                else
                {
                    return FFEA_ERROR;
                }
            }

            force_mag = spring_array[i].k * (mag<scalar, arr3>(sep.data) - spring_array[i].l);
            force0.x = force_mag * sep_norm.x;
            force0.y = force_mag * sep_norm.y;
            force0.z = force_mag * sep_norm.z;

            force1.x = -1 * force_mag * sep_norm.x;
            force1.y = -1 * force_mag * sep_norm.y;
            force1.z = -1 * force_mag * sep_norm.z;
            active_blob_array[spring_array[i].blob_index[0]]->add_force_to_node(force0, spring_array[i].node_index[0]);
            active_blob_array[spring_array[i].blob_index[1]]->add_force_to_node(force1, spring_array[i].node_index[1]);
        }
    }
    return FFEA_OK;
}

scalar World::get_spring_field_energy(int index0, int index1)
{

    // Sum over all field
    if (index0 == -1 || index1 == -1)
    {
        scalar energy = 0.0;
        for (int i = 0; i < params.num_blobs; ++i)
        {
            for (int j = 0; j < params.num_blobs; ++j)
            {
                energy += springfieldenergy[i][j];
            }
        }

        return energy;
    }
    else if (index0 == index1)
    {
        return springfieldenergy[index0][index1];
    }
    else
    {

        // Order of blob indices is unknown in the calculations, so must add
        return springfieldenergy[index0][index1] + springfieldenergy[index1][index0];
    }
}

int World::get_next_script_tag(FILE *in, char *buf)
{
    if (fscanf(in, "%*[^<]") != 0)
    {
        printf("White space removal error in get_next_script_tag(). Something odd has happened: this error should never occur...\n");
    }
    if (fscanf(in, "<%255[^>]>", buf) != 1)
    {
        FFEA_error_text();
        printf("Error reading tag in script file.\n");
        return FFEA_ERROR;
    }

    return FFEA_OK;
}

void World::apply_dense_matrix(scalar *y, scalar *M, scalar *x, int N)
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        y[i] = 0;
        for (j = 0; j < N; j++)
            y[i] += M[i * N + j] * x[j];
    }
}

void World::do_es()
{
    printf("Building BEM matrices\n");
    PB_solver.build_BEM_matrices();

    //	PB_solver.print_matrices();

    // Build the poisson matrices for each blob
    printf("Building Poisson matrices\n");
    for (int i = 0; i < params.num_blobs; i++)
    {
        active_blob_array[i]->build_poisson_matrices();
    }

    printf("Solving\n");
    for (int dual = 0; dual < 30; dual++)
    {

        // Perform Poisson solve step (K phi + E J_Gamma = rho)
        // Obtain the resulting J_Gamma
        int master_index = 0;
        for (int i = 0; i < params.num_blobs; i++)
        {
            active_blob_array[i]->solve_poisson(&phi_Gamma[master_index], &J_Gamma[master_index]);
            master_index += active_blob_array[i]->get_num_faces();
        }

        // Apply -D matrix to J_Gamma vector (work_vec = -D J_gamma)
        PB_solver.get_D()->apply(J_Gamma, work_vec);
        for (int i = 0; i < total_num_surface_faces; i++)
        {
            work_vec[i] *= -1;
        }

        // Solve for C matrix (C phi_Gamma = work_vec)
        nonsymmetric_solver.solve(PB_solver.get_C(), phi_Gamma, work_vec);

        scalar sum = 0.0, sumj = 0.0;
        for (int i = 0; i < total_num_surface_faces; i++)
        {
            sum += phi_Gamma[i];
            sumj += J_Gamma[i];
        }
        printf("\n");

        printf("<yophi> = %e <J> = %e\n", sum / total_num_surface_faces, sumj / total_num_surface_faces);

        //		for(i = 0; i < total_num_surface_faces; i++) {
        //			printf("<JPBSN> %d %e\n", dual, J_Gamma[i]);
        //		}
    }

    //	scalar sum = 0.0, sumj = 0.0;
    //	for(i = 0; i < total_num_surface_faces; i++) {
    //		sum += phi_Gamma[i];
    //		sumj += J_Gamma[i];
    //	}
    //	printf("<phi> = %e <J> = %e\n", sum/total_num_surface_faces, sumj/total_num_surface_faces);
    //
    //	for(i = 0; i < total_num_surface_faces; i++) {
    //		printf("<J> %e\n", J_Gamma[i]);
    //	}

    //	blob_array[0].print_phi();
}

/**
 * @brief Writes eigensystems to files in order
 * @details This function writes the eigenvalues calculated for enms and dmms to files
 * */

void World::write_eig_to_files(scalar *evals_ordered, scalar **evecs_ordered, int num_modes, int num_nodes)
{

    // Get some filenames
    vector<string> all;
    string val_out_fname, vec_out_fname, base, ext;
    boost::split(all, params.trajectory_out_fname, boost::is_any_of("."));
    ext = "." + all.at(all.size() - 1);
    base = boost::erase_last_copy(params.trajectory_out_fname, ext);
    val_out_fname = base + ".evals";
    vec_out_fname = base + ".evecs";

    // Open the files and write the modes out
    FILE *valfout, *vecfout;
    valfout = fopen(val_out_fname.c_str(), "wb");
    vecfout = fopen(vec_out_fname.c_str(), "wb");
    for (int i = 0; i < 3 * num_nodes; ++i)
    {
        if (i < num_modes)
        {
            fprintf(valfout, "%6.3e\n", evals_ordered[i]);
        }
        for (int j = 0; j < num_modes; ++j)
        {
            fprintf(vecfout, "%6.3e ", evecs_ordered[j][i]);
        }
        fprintf(vecfout, "\n");
    }
    fclose(valfout);
    fclose(vecfout);
}

/**
 * @brief Calculates a pseudo-trajectory by varying an eigenvector by a constant factor
 * @details This function takes an Eigen::VectorXd and applies it as a series of translations
 * to the given blob.
 * */

void World::make_trajectory_from_eigenvector(string traj_out_fname, int blob_index, int mode_index, Eigen_VectorX evec, scalar step)
{

    int i, j, from_index = 0, to_index = 0;
    scalar dx;

    // Convert weird eigen thing into a nice vector3 list
    vector<vector3> node_trans(active_blob_array[blob_index]->get_num_linear_nodes());

    // Open file
    FILE *fout;
    fout = fopen(traj_out_fname.c_str(), "wb");

    // Header Stuff
    fprintf(fout, "FFEA_trajectory_file\n\nInitialisation:\nNumber of Blobs 1\nNumber of Conformations 1\nBlob 0:	Conformation 0 Nodes %d\n\n*\n", active_blob_array[blob_index]->get_num_nodes());

    // Initial centroid, to move around
    //active_blob_array[blob_index]->position(0,0,0);
    for (i = 0; i < 21; ++i)
    {

        /* Build a frame */

        // Get eigenvector multiplier
        if (i == 0)
        {
            dx = 0.0;
        }
        else if (i > 0 && i < 5)
        {
            dx = step;
        }
        else if (i > 5 && i <= 15)
        {
            dx = -step;
        }
        else
        {
            dx = step;
        }

        // Get some node translations
        int num_linear_nodes = active_blob_array[blob_index]->get_num_linear_nodes();
        for (j = 0; j < num_linear_nodes; ++j)
        {
            node_trans[j].x = evec[3 * j + 0] * dx;
            node_trans[j].y = evec[3 * j + 1] * dx;
            node_trans[j].z = evec[3 * j + 2] * dx;
        }

        // Translate all the nodes
        active_blob_array[blob_index]->translate_linear(node_trans.data());
        fprintf(fout, "Blob 0, Conformation 0, step %d\n", i);
        active_blob_array[blob_index]->write_nodes_to_file(fout);
        fprintf(fout, "*\n");
        print_trajectory_conformation_changes(fout, i, &from_index, &to_index);
    }
    fclose(fout);
}

void World::print_evecs_to_file(string fname, Eigen_MatrixX ev, int num_rows, int num_modes)
{

    int i, j;
    FILE *fout;
    fout = fopen(fname.c_str(), "wb");

    // Skip the zero modes
    for (i = 6; i < num_modes + 6; ++i)
    {
        for (j = 0; j < num_rows; ++j)
        {
            fprintf(fout, "%6.3f ", ev.col(i)[j]);
        }
        fprintf(fout, "\n");
    }
    //for(i = num_modes + 5; i > 5; --i) {
    //    for(j = 0; j < num_rows; ++j) {
    //        fprintf(fout, "%6.3f ", ev.col(i)[j]);
    //    }
    //    fprintf(fout, "\n");
    //}
    fclose(fout);
}

void World::print_evals_to_file(string fname, Eigen_VectorX ev, int num_modes, scalar scale)
{

    int i;
    FILE *fout;
    fout = fopen(fname.c_str(), "wb");

    // Skip the zero modes
    for (i = 6; i < num_modes + 6; ++i)
    {
        fprintf(fout, "%6.3e\n", ev[i] * scale);
    }
    //for(i = num_modes + 5; i > 5; --i) {
    //    fprintf(fout, "%6.3e\n", ev[i] * scale);
    //}
    fclose(fout);
}

void World::write_output_header(FILE *fout, string fname)
{

    // Write all header data. script fname, time and date etc
    fprintf(fout, "FFEA Global Measurement File\n\nSimulation Details:\n");

    print_ffea_version(measurement_out);
    print_ffea_compilation_details(measurement_out);

    time_t now = time(0);
    tm *ltm = localtime(&now);
    fprintf(fout, "\tSimulation Began on %d/%d/%d at %d:%d:%d\n", ltm->tm_mday, 1 + ltm->tm_mon, 1900 + ltm->tm_year, ltm->tm_hour, ltm->tm_min, ltm->tm_sec);
    fprintf(fout, "\tScript Filename = %s\n", fname.c_str());
    fprintf(fout, "\tSimulation Type = %s\n\n", "Full");
}

void World::write_pre_print_to_trajfile(int step)
{
    for (int i = 0; i < params.num_blobs; i++)
    {
        fprintf(trajectory_out, "Blob %d, Conformation %d, step %d\n", i, active_blob_array[i]->get_conformation_index(), step);
        active_blob_array[i]->write_pre_print_to_file(trajectory_out);
    }
    // Mark completed end of step with an asterisk (so that the restart code will know if this is a fully written step or if it was cut off half way through due to interrupt)
    fprintf(trajectory_out, "*\n");

    // And print the states.
    fprintf(trajectory_out, "Conformation Changes:\n");
    for (int i = 0; i < params.num_blobs; ++i)
    {
        if (params.calc_kinetics == 1 && active_blob_array[i]->toBePrinted_state[0] != active_blob_array[i]->toBePrinted_state[1])
        {
            printf("\tBlob %d - Conformation %d -> Conformation %d\n", i, active_blob_array[i]->toBePrinted_conf[0], active_blob_array[i]->toBePrinted_conf[1]);
            printf("\t		State %d -> State %d\n", active_blob_array[i]->toBePrinted_state[0], active_blob_array[i]->toBePrinted_state[1]);
        }

        // Print to file
        fprintf(trajectory_out, "Blob %d: Conformation %d -> Conformation %d\n", i, active_blob_array[i]->toBePrinted_conf[0], active_blob_array[i]->toBePrinted_conf[1]);
    }
    fprintf(trajectory_out, "*\n");

    fflush(trajectory_out);
}

/** Write trajectory for each blob, then do blob specific measurements (which are needed for globals, but only explicitly printed if "-d" was used) */
void World::print_trajectory_and_measurement_files(int step, scalar wtime)
{

    // ONSCREEN progress:
    if (step % (params.check * 10) != 0)
    {
        printf("\rstep = %d", step);
        fflush(stdout);
    }
    else
    {
        printf("\rstep = %d (simulation time = %.2fns, wall clock time = %.3f hrs)\n", step, step * params.dt * (mesoDimensions::time / 1e-9), (omp_get_wtime() - wtime) / 3600.0);
    }

    // TRAJECTORY file: can be printed serially, or in parallel:
#ifdef FFEA_PARALLEL_FUTURE
    // TRAJECTORY PARALLEL:
    // CHECK // CHECK // CHECK //
    bool its = false;
    if (thread_writingTraj.valid())
    {
        if (thread_writingTraj.wait_for(std::chrono::milliseconds(0)) == future_status::ready)
        {
            its = true;
        }
    }
    if (!its)
        cout << " We're waiting while writing the trajectory" << endl;
    // END CHECK // END CHECK // END CHECK //
    thread_writingTraj.get();
#ifdef FFEA_PARALLEL_PER_BLOB
#pragma omp parallel for default(none) shared(step) schedule(static)
#endif
    for (int i = 0; i < params.num_blobs; i++)
    {
        // store the node data for this blob
        active_blob_array[i]->pre_print();
    }
    thread_writingTraj = std::async(std::launch::async, &World::write_pre_print_to_trajfile, this, step);
    // write_pre_print_to_trajfile(step); // serial version
#else
    // TRAJECTORY SERIAL:
    for (int i = 0; i < params.num_blobs; i++)
    {

        // Write the node data for this blob
        fprintf(trajectory_out, "Blob %d, Conformation %d, step %d\n", i, active_blob_array[i]->get_conformation_index(), step);
        active_blob_array[i]->write_nodes_to_file(trajectory_out);
    }
    // Mark completed end of step with an asterisk (so that the restart code will know if this is a fully written step or if it was cut off half way through due to interrupt)
    fprintf(trajectory_out, "*\n");

    // And print the states.
    fprintf(trajectory_out, "Conformation Changes:\n");
    for (int i = 0; i < params.num_blobs; ++i)
    {
        if (params.calc_kinetics == 1 && active_blob_array[i]->get_previous_state_index() != active_blob_array[i]->get_state_index())
        {
            printf("\tBlob %d - Conformation %d -> Conformation %d\n", i, active_blob_array[i]->get_previous_conformation_index(), active_blob_array[i]->get_conformation_index());
            printf("\t		State %d -> State %d\n", active_blob_array[i]->get_previous_state_index(), active_blob_array[i]->get_state_index());
        }

        // Print to file
        fprintf(trajectory_out, "Blob %d: Conformation %d -> Conformation %d\n", i, active_blob_array[i]->get_previous_conformation_index(), active_blob_array[i]->get_conformation_index());
    }
    fprintf(trajectory_out, "*\n");

    // Force print in case of ctrl + c stop
    fflush(trajectory_out);
#endif
    // TRAJECTORY END

    //Write rod trajectory (skip first frame if this is a restart)
    for (int i = 0; i < params.num_rods; i++)
    {
        if (rod_array[i]->restarting)
        {
            rod_array[i]->restarting = false;
        }
        else
        {
            rod_array[i]->write_frame_to_file();
        }
    }

    // Detailed Measurement Stuff.
    // Stuff needed on each blob, and in global energy files
    if (detailed_meas_out != NULL)
    {
        fprintf(detailed_meas_out, "%-14.6e", step * params.dt * mesoDimensions::time);
    }

#ifdef FFEA_PARALLEL_PER_BLOB
#pragma omp parallel for default(none) schedule(static)
#endif
    for (int i = 0; i < params.num_blobs; i++)
    {
        // Calculate properties for this blob
        active_blob_array[i]->make_measurements();
    }

    // If necessary, write this stuff to a separate file
    for (int i = 0; i < params.num_blobs; i++)
    {
        if (detailed_meas_out != NULL)
        {
            active_blob_array[i]->write_measurements_to_file(detailed_meas_out);
        }
    }

    // Global Measurement Stuff
    make_measurements();
    write_measurements_to_file(measurement_out, step);

    if (detailed_meas_out != NULL)
    {
        write_detailed_measurements_to_file(detailed_meas_out);
    }
    fflush(measurement_out);

    if (params.trajbeads_fname_set == 1)
    {
        pc_solver.write_beads_to_file(trajbeads_out, step);
    }
}

/** Print the RNG values responsible for the current state of the system i.e. do this before all forces are calculated, which advances the system */
void World::print_checkpoints()
{

    // CHECKPOINT - Write the state of the RNGs:
    // REWIND!
    rewind(checkpoint_out);
    // Header for the thermal stresses:
    int thermal_seeds = num_seeds;
    if (params.calc_kinetics == 1)
        thermal_seeds -= 1;
    fprintf(checkpoint_out, "RNGStreams dedicated to the thermal stress: %d\n", thermal_seeds);
    //cout << "hi" << endl << flush;
    unsigned long state[6];
    // First save the state of the running threads:
    for (int i = 0; i < num_threads; i++)
    {
        rng[i].GetState(state);
        fprintf(checkpoint_out, "%lu %lu %lu %lu %lu %lu\n", state[0], state[1], state[2],
                state[3], state[4], state[5]);
        //for(int j = 0; j < 6; ++j) {
        // cout << " " << state[j];
        //}
        // cout << endl;
    }

    // If there were more threads running on the previous run, we'll save them too:
    int oldThreads = thermal_seeds - num_threads;
    //cout << oldThreads << " " << thermal_seeds << " " << num_threads << endl << flush;
    for (int i = 0; i < oldThreads; i++)
    {
        fprintf(checkpoint_out, "%lu %lu %lu %lu %lu %lu\n", Seeds[i + num_threads][0],
                Seeds[i + num_threads][1], Seeds[i + num_threads][2], Seeds[i + num_threads][3],
                Seeds[i + num_threads][4], Seeds[i + num_threads][5]);
        // for(int j = 0; j < 6; ++j) {
        //    cout << " " << Seeds[i+num_threads][j];
        //     }
        //  cout << endl;
    }
    //cout << "hi" << endl << flush;
    // If we're doing kinetics, we're saving the state of the extra RNG:
    if (params.calc_kinetics)
    {
        fprintf(checkpoint_out, "RNGStream dedicated to the kinetics:\n");
        kinetic_rng->GetState(state);
        fprintf(checkpoint_out, "%lu %lu %lu %lu %lu %lu\n", state[0], state[1], state[2],
                state[3], state[4], state[5]);
    }
    //cout << "hi" << endl << flush;
    fflush(checkpoint_out);
    // Done with the checkpoint!
}

void World::make_measurements()
{

    int i, j, total_num_nodes = 0;

    // Set stuff to zero
    kineticenergy = 0.0;
    strainenergy = 0.0;
    springenergy = 0.0;
    ssintenergy = 0.0;
    preCompenergy = 0.0;
    rmsd = 0.0;
    vector3_set_zero(CoG);

    vector3 bCoG;

    // Sum stuff from blobs
    for (i = 0; i < params.num_blobs; ++i)
    {
        kineticenergy += active_blob_array[i]->get_kinetic_energy();
        strainenergy += active_blob_array[i]->get_strain_energy();
        rmsd += (active_blob_array[i]->get_rmsd() * active_blob_array[i]->get_rmsd()) * active_blob_array[i]->get_num_nodes();
        active_blob_array[i]->get_stored_centroid(bCoG.data);
        CoG.x += bCoG.x * active_blob_array[i]->get_num_nodes();
        CoG.y += bCoG.y * active_blob_array[i]->get_num_nodes();
        CoG.z += bCoG.z * active_blob_array[i]->get_num_nodes();

        total_num_nodes += active_blob_array[i]->get_num_nodes();
    }

    // And divide by num_nodes
    rmsd = sqrt(rmsd / total_num_nodes);
    arr3Resize<scalar, arr3>(1.0 / total_num_nodes, CoG.data);

    // Now global stuff
    vector3 a, b, c;
    if (params.calc_springs != 0)
    {
        for (i = 0; i < params.num_blobs; ++i)
        {
            for (j = 0; j < params.num_blobs; ++j)
            {
                springfieldenergy[i][j] = 0.0;
            }
        }

        for (i = 0; i < num_springs; ++i)
        {
            active_blob_array[spring_array[i].blob_index[0]]->get_node(spring_array[i].node_index[0], a.data);
            active_blob_array[spring_array[i].blob_index[1]]->get_node(spring_array[i].node_index[1], b.data);
            arr3arr3Substract<scalar, arr3>(a.data, b.data, c.data);
            springfieldenergy[spring_array[i].blob_index[0]][spring_array[i].blob_index[1]] += 0.5 * spring_array[i].k * (mag<scalar, arr3>(c.data) - spring_array[i].l) * (mag<scalar, arr3>(c.data) - spring_array[i].l);
        }

        springenergy = get_spring_field_energy(-1, -1);
    }

    if (params.calc_ssint == 1 || params.calc_steric == 1)
    {
        ssintenergy = vdw_solver->get_field_energy(-1, -1);
    }

    if (params.calc_preComp == 1)
    {
        preCompenergy = pc_solver.get_field_energy(-1, -1);
    }
}

void World::write_measurements_to_file(FILE *fout, int step)
{

    // In same order as initialisation
    fprintf(fout, "%-14.6e", step * params.dt * mesoDimensions::time);
    if (mass_in_system)
    {
        fprintf(fout, "%-14.6e", kineticenergy * mesoDimensions::Energy);
    }
    fprintf(fout, "%-14.6e", strainenergy * mesoDimensions::Energy);
    fprintf(fout, "%-14.6e%-14.6e%-14.6e", CoG.x * mesoDimensions::length, CoG.y * mesoDimensions::length, CoG.z * mesoDimensions::length);
    fprintf(fout, "%-14.6e", rmsd * mesoDimensions::length);
    if (params.calc_springs != 0)
    {
        fprintf(fout, "%-14.6e", springenergy * mesoDimensions::Energy);
    }
    if (params.calc_ssint == 1 || params.calc_steric == 1)
    {
        fprintf(fout, "%-15.6e", ssintenergy * mesoDimensions::Energy);
    }
    if (params.calc_preComp != 0)
    {
        fprintf(fout, "%-14.6e", preCompenergy * mesoDimensions::Energy);
    }

    fprintf(fout, "\n");
    fflush(fout);
}

void World::write_detailed_measurements_to_file(FILE *fout)
{

    // In same order as initialisation
    int i, j;
    for (i = 0; i < params.num_blobs; ++i)
    {
        for (j = i; j < params.num_blobs; ++j)
        {
            // White space for blob index bit
            fprintf(fout, "       ");
            if (active_blob_array[i]->there_is_ssint() && active_blob_array[j]->there_is_ssint())
            {
                fprintf(detailed_meas_out, "%-15.6e", vdw_solver->get_field_energy(i, j) * mesoDimensions::Energy);
            }
            if (active_blob_array[i]->there_are_springs() && active_blob_array[j]->there_are_springs() * mesoDimensions::Energy)
            {
                fprintf(detailed_meas_out, "%-14.6e", get_spring_field_energy(i, j) * mesoDimensions::Energy);
            }
            if (active_blob_array[i]->there_are_beads() && active_blob_array[j]->there_are_beads() * mesoDimensions::Energy)
            {
                fprintf(detailed_meas_out, "%-14.6e", pc_solver.get_field_energy(i, j) * mesoDimensions::Energy);
            }
        }
    }
    fprintf(fout, "\n");
    fflush(fout);
}

void World::print_trajectory_conformation_changes(FILE *fout, int step, int *from_index, int *to_index)
{

    // Check input
    int *to;
    int *from;
    if (to_index == NULL || from_index == NULL)
    {
        to = new int[params.num_blobs];
        from = new int[params.num_blobs];
        for (int i = 0; i < params.num_blobs; ++i)
        {
            to[i] = active_blob_array[i]->conformation_index;
            from[i] = active_blob_array[i]->conformation_index;
        }
    }
    else
    {
        to = to_index;
        from = from_index;
    }

    // Inform whoever is watching of changes
    if (params.calc_kinetics == 1 && (step - 1) % params.kinetics_update == 0)
    {
        printf("Conformation Changes:\n");
    }
    fprintf(fout, "Conformation Changes:\n");
    for (int i = 0; i < params.num_blobs; ++i)
    {
        if (params.calc_kinetics == 1 && (step - 1) % params.kinetics_update == 0)
        {
            printf("\tBlob %d - Conformation %d -> Conformation %d\n", i, from[i], to[i]);
        }

        // Print to file
        fprintf(fout, "Blob %d: Conformation %d -> Conformation %d\n", i, from[i], to[i]);
    }
    fprintf(fout, "*\n");

    // Force print in case of ctrl + c stop
    fflush(fout);

    // Free data
    if (to_index == NULL || from_index == NULL)
    {
        delete[] to;
        delete[] from;
    }
}

void World::print_kinetic_files(int step)
{

    if (params.calc_kinetics == 0)
        return;

    // Print to specific file
    if (kinetics_out != NULL)
    {
        fprintf(kinetics_out, "%d", step);
        for (int i = 0; i < params.num_blobs; ++i)
        {
            fprintf(kinetics_out, " %d %d", active_blob_array[i]->get_state_index(), active_blob_array[i]->get_conformation_index());
        }
        fprintf(kinetics_out, "\n");
        fflush(kinetics_out);
    }

    // And now previous state is the current state
    for (int i = 0; i < params.num_blobs; ++i)
    {
        active_blob_array[i]->set_previous_state_index(active_blob_array[i]->get_state_index());
        active_blob_array[i]->set_previous_conformation_index(active_blob_array[i]->get_conformation_index());
    }
}

void World::print_static_trajectory(int step, scalar wtime, int blob_index)
{
    printf("Printing single trajectory of Blob %d for viewer\n", blob_index);
    // Write the node data for this blob
    fprintf(trajectory_out, "Blob %d, step %d\n", blob_index, step);
    active_blob_array[blob_index]->write_nodes_to_file(trajectory_out);
}

void World::calc_blob_corr_matrix(int num_blobs, scalar *blob_corr)
{

    //calculates blob corrections for periodic interactions
    vector3 com, com2;
    for (int i = 0; i < num_blobs; ++i)
    {
        //sets blob distance from itself to 0
        blob_corr[i * num_blobs * 3 + i * 3] = 0;
        blob_corr[i * num_blobs * 3 + i * 3 + 1] = 0;
        blob_corr[i * num_blobs * 3 + i * 3 + 2] = 0;

        active_blob_array[i]->get_stored_centroid(com.data);
        for (int j = i + 1; j < num_blobs; ++j)
        {

            active_blob_array[j]->get_stored_centroid(com2.data);
            blob_corr[i * num_blobs * 3 + j * 3] = box_dim.x * floor((com2.x - com.x + 0.5 * box_dim.x) / box_dim.x);
            blob_corr[i * num_blobs * 3 + j * 3 + 1] = box_dim.y * floor((com2.y - com.y + 0.5 * box_dim.y) / box_dim.y);
            blob_corr[i * num_blobs * 3 + j * 3 + 2] = box_dim.z * floor((com2.z - com.z + 0.5 * box_dim.z) / box_dim.z);

            blob_corr[j * num_blobs * 3 + i * 3] = (-1) * blob_corr[i * num_blobs * 3 + j * 3];
            blob_corr[j * num_blobs * 3 + i * 3 + 1] = (-1) * blob_corr[i * num_blobs * 3 + j * 3 + 1];
            blob_corr[j * num_blobs * 3 + i * 3 + 2] = (-1) * blob_corr[i * num_blobs * 3 + j * 3 + 2];
        }
    }
}

void World::do_nothing()
{
    // that means nothing.
}

int World::die_with_dignity(int step, scalar wtime)
{

    FFEA_error_text();
    printf("A problem occurred when...\n");
    printf("Simulation ran for %2f seconds (wall clock time) before error ocurred\n", (omp_get_wtime() - wtime));
    printf("Detected fatal errors in this system update. Exiting now...\n");

    // attempt to print out the final (bad) time step (if step != step_initial)
    if (step != step_initial)
    {
        printf("Dumping final step:\n");
        print_trajectory_and_measurement_files(step, wtime);
        print_kinetic_files(step);
    }

    return FFEA_ERROR;
}

#ifdef FFEA_PARALLEL_FUTURE
int World::prebuild_nearest_neighbour_lookup_wrapper(scalar cell_size)
{
    return lookup.prebuild_nearest_neighbour_lookup(cell_size);
}

bool World::updatingVdWLL()
{
    return thread_updatingVdWLL.valid();
}

bool World::updatingPCLL()
{
    return thread_updatingPCLL.valid();
}

bool World::updatingVdWLL_ready_to_swap()
{

    bool its = false;
    if (thread_updatingVdWLL.valid())
    {
        if (thread_updatingVdWLL.wait_for(std::chrono::milliseconds(0)) == future_status::ready)
        {
            its = true;
        }
    }
    return its;
}

bool World::updatingPCLL_ready_to_swap()
{

    bool its = false;
    if (thread_updatingPCLL.valid())
    {
        if (thread_updatingPCLL.wait_for(std::chrono::milliseconds(0)) == future_status::ready)
        {
            its = true;
        }
    }
    return its;
}

int World::catch_thread_updatingVdWLL(int step, scalar wtime, int where)
{

    if (updatingVdWLL() == false)
    { // i. e., thread has been already catched!
        cout << "trying to catch from: " << where << ", but updatingVdWLL was false" << endl;
        return 0;
    }

    if (thread_updatingVdWLL.get() == FFEA_ERROR)
    {
        return die_with_dignity(step, wtime);
    }
    if (lookup.safely_swap_layers() == FFEA_ERROR)
    {
        return die_with_dignity(step, wtime);
    }

    return FFEA_OK;
}

int World::catch_thread_updatingPCLL(int step, scalar wtime, int where)
{

    if (updatingPCLL() == false)
    { // i. e., thread has been already catched!
        cout << "trying to catch from: " << where << ", but updatingPCLL was false" << endl;
        return 0;
    }

    if (thread_updatingPCLL.get() == FFEA_ERROR)
    { // we'll need to change this check
        return die_with_dignity(step, wtime);
    }
    if (pc_solver.safely_swap_pc_layers() == FFEA_ERROR)
    {
        return die_with_dignity(step, wtime);
    }

    return FFEA_OK;
}

#endif

// Well done for reading this far! Hope this makes you smile.

/*
quu..__
 $$$b  `---.__
  "$$b        `--.                          ___.---uuudP
   `$$b           `.__.------.__     __.---'      $$$$"              .
     "$b          -'            `-.-'            $$$"              .'|
       ".                                       d$"             _.'  |
         `.   /                              ..."             .'     |
           `./                           ..::-'            _.'       |
            /                         .:::-'            .-'         .'
           :                          ::''\          _.'            |
          .' .-.             .-.           `.      .'               |
          : /'$$|           .@"$\           `.   .'              _.-'
         .'|$u$$|          |$$,$$|           |  <            _.-'
         | `:$$:'          :$$$$$:           `.  `.       .-'
         :                  `"--'             |    `-.     \
        :##.       ==             .###.       `.      `.    `\
        |##:                      :###:        |        >     >
        |#'     `..'`..'          `###'        x:      /     /
         \                                   xXX|     /    ./
          \                                xXXX'|    /   ./
          /`-.                                  `.  /   /
         :    `-  ...........,                   | /  .'
         |         ``:::::::'       .            |<    `.
         |             ```          |           x| \ `.:``.
         |                         .'    /'   xXX|  `:`M`M':.
         |    |                    ;    /:' xXXX'|  -'MMMMM:'
         `.  .'                   :    /:'       |-'MMMM.-'
          |  |                   .'   /'        .'MMM.-'
          `'`'                   :  ,'          |MMM<
            |                     `'            |tbap\
             \                                  :MM.-'
              \                 |              .''
               \.               `.            /
                /     .:::::::.. :           /
               |     .:::::::::::`.         /
               |   .:::------------\       /
              /   .''               >::'  /
              `',:                 :    .'
*/
