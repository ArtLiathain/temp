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

#ifndef PRECOMP_H_INCLUDED
#define PRECOMP_H_INCLUDED

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <algorithm>

#include "FFEA_return_codes.h"
#include "dimensions.h"
#include "FFEA_user_info.h"
#include "mat_vec_types.h"
#include "LinkedListCube.h"
#include "SimulationParams.h"

// WARNING: Blob.h will be included after defining PreComp_params! 

using namespace std;

#include "Blob.h"

class PreComp_solver{
public:
  PreComp_solver();
  ~PreComp_solver();
  int init(PreComp_params *pc_params, SimulationParams *params, Blob **blob_array);
  int solve(scalar *blob_corr=NULL); ///< calculate the forces using a straightforward double loop.
  int solve_using_neighbours();  ///< calculate the forces using linkedlists.
  int solve_using_neighbours_non_critical(scalar *blob_corr=NULL);  ///< using linkedlists, calculate twice the forces to avoid any critical regions.
  void reset_fieldenergy(); 
  scalar get_U(scalar x, int typei, int typej);
  scalar get_F(scalar x, int typei, int typej);
  scalar get_field_energy(int index0, int index1);

  int compute_bead_positions(); ///< calculate b_pos, the absolute positions of the beads. 

  int build_pc_nearest_neighbour_lookup(); ///< put the beads on the grid.
  int prebuild_pc_nearest_neighbour_lookup_and_swap(); ///< put the beads on the grid.
  int prebuild_pc_nearest_neighbour_lookup(); ///< put the beads on the grid.
  int safely_swap_pc_layers(); ///< swap the two LinkedLists. 

  void write_beads_to_file(FILE *fout, int timestep); ///< write beads to file, for the current timestep

private: 
  /** msgc and msg are helpful while developing */
  int msgc;
  int msg(string whatever); 
  int msg(int whatever); 
  
  int read_tabulated_values(PreComp_params &pc_params, string kind, scalar *Z, scalar scale_Z);

  int calc_force_from_pot();
  
  scalar finterpolate(scalar *Z, scalar x, int typei, int typej);

  // stuff related to the LinkedLists:
  LinkedListCube<int> pcLookUp; ///< the linkedlist itself
  scalar pcVoxelSize;    ///< the size of the voxels.
  int pcVoxelsInBox[3];  ///< num of voxels per side.
  static const int adjacent_cells[27][3]; 
  

  /** delta x in tabulated potentials and forces"  */
  scalar Dx; 
  /** x_range */ 
  scalar x_range[2];
  /** squared x_range */ 
  scalar x_range2[2];
  /** number of pre-computed values per table */ 
  int n_values; 
  /** total number of type interactions */
  int nint; 
  /** number of types of "beads */
  int ntypes; 
  /** pointer to array containing all the values for all the pair potentials. */
  scalar *U;
  /** pointer to array containing all the values for all the pair forces. */
  scalar *F;
  /** interacting elements */
  typedef tetra_element_linear* TELPtr;
  TELPtr *b_elems;
  TELPtr *b_unq_elems; 
  /** bead types */
  int *b_types; 
  /** map b_unq_elems to beads */
  int *map_e_to_b;
  /** number of beads */ 
  int n_beads; 
  /** number of different elements */
  int num_diff_elems; 
  /** relative position of the beads to the element they belong, xyzxyzxyz... */
  scalar *b_rel_pos;
  /** absolute position of the beads */
  scalar *b_pos; 
  /** forces to be applied */
  scalar *b_forces;
  /** list of the daddy blob */
  int *b_daddyblob; 
  /** bool "matrix" (array) storing for every pair if it is active or not */
  bool *isPairActive;

  /** variables stypes, b_elems_ndx, and b_blob_ndx will only be used if writing traj */ 
  vector<string> stypes; ///< string types for the beads; sorry it is a c++ vector
  int *b_elems_ndx; ///< array with the corresponding element index.
  int *b_blob_ndx; ///< array with the corresponding blob index.
  

  /** Variables to store the energy field data between each pair of blobs */
  scalar *fieldenergy;
  int num_blobs;
  int num_threads;

};

#endif
