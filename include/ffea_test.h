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

/*
 *      rod_blob_interface.cpp
 *	Author: Rob Welch, University of Leeds
 *	Email: py12rw@leeds.ac.uk

 *	Author: Ryan Cocking, University of Leeds
 *	Email: bsrctb@leeds.ac.uk
 */

#include <iostream>
#include <omp.h>
#include <cstdlib>

#include "FFEA_return_codes.h"
#include "FFEA_user_info.h"
#include "mat_vec_types.h"
#include "mesh_node.h"
#include "Blob.h"
#include "Face.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseSubstitutionSolver.h"
#include "World.h"
#include "rod_interactions.h"


// ffea_test allows you to build unit tests in FFEA that can be accessed via a
// .ffeatest script file The reason is that otherwise, you'd have to compile a
// new binary, and you can't easily separate out a small section of FFEA
// codebase to compile - you basically have to compile the whole binary each
// time.

struct ffea_test
{
    static int do_ffea_test(std::string filename);
    static int connection_test();
    static int arbitrary_equilibrium_twist();
    static int connection_orientation_test();
    static int arbitrary_equilibrium_bend();
    static int identify_face();
    static int connection_energy();
    static int connection_energy_2();
    static int jacobian_rotate();
    static int connection_energy_3();
    static int connection_propagation_every_way();
    static int connection_propagation(int mode, bool ends_at_rod);
    static int recover_normal();
    static int dump_twist_info();
    static int lower_sphere();
    static int point_lies_within_rod_element();
    static int line_connecting_rod_elements();
    static int rod_neighbour_list_construction();
    static int rod_steric_lj_potential();
    static int nearest_image_pbc();
    static int rod_vdw_site_placement();
};
