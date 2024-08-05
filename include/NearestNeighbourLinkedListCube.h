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

#ifndef NEARESTNEIGHBOURLINKEDLISTCUBE_H_INCLUDED
#define NEARESTNEIGHBOURLINKEDLISTCUBE_H_INCLUDED

#include <cmath>
#include "mat_vec_types.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "LinkedListCube.h"
#include "Face.h"

class NearestNeighbourLinkedListCube : public LinkedListCube<Face> {
public:
    /** Build the nearest neighbour look up cube given the spatial cell size */
    int build_nearest_neighbour_lookup(scalar h);

    /** Build the nearest neighbour look up cube given the spatial cell size */
    int prebuild_nearest_neighbour_lookup_and_swap(scalar h);

    /** Build the nearest neighbour look up cube given the spatial cell size */
    int prebuild_nearest_neighbour_lookup(scalar h);
};

#endif
