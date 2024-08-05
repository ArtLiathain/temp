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

#ifndef GAUSSIANQUADRATURE_1D_H_INCLUDED
#define GAUSSIANQUADRATURE_1D_H_INCLUDED

#include <cmath>
#include "mat_vec_types.h"

typedef struct {
    /** Weight */
    scalar W;

    /** Abscissa */
    scalar zeta;
} gauss_point_1d;


static gauss_point_1d gq_line[] = {
    /* 1 point - precision 1 */
    {2.0, 0.0},

    /* 2 point - precision 2 */
    {1.0, 1.0 / sqrt(3)},
    {1.0, -1.0 / sqrt(3)},

    /* 3 point - precision 3 */
    {8.0 / 9.0, 0.0},
    {5.0 / 9.0, sqrt(3.0 / 5.0)},
    {5.0 / 9.0, -sqrt(3.0 / 5.0)},

    /* 4 point - precision 4 */
    {0.652145154862546, 0.339981043584856},
    {0.652145154862546, -0.339981043584856},
    {0.347854845137454, 0.861136311594053},
    {0.347854845137454, -0.861136311594053},

    /* 5 point - precision 5 */
    {0.568888888888889, 0.0},
    {0.478628670499366, 0.538469310105683},
    {0.478628670499366, -0.538469310105683},
    {0.236926885056189, 0.906179845938664},
    {0.236926885056189, -0.906179845938664},

    /* 6 Point - precision 6 */
    {0.467913934572691, 0.238619186083197},
    {0.467913934572691, -0.238619186083197},
    {0.360761573048139, 0.661209386466265},
    {0.360761573048139, -0.661209386466265},
    {0.171324492379170, 0.932469514203152},
    {0.171324492379170, -0.932469514203152}

};

static int gq_precision_index[] = {
    0, // Ignore this line (no precision 0)
    0, // Precision 1: Starts at index 0, 1 point
    1, // Precision 2: Starts at index 1, 2 point
    3, // Precision 3: Starts at index 3, 3 points
    6, // Precision 4: Starts at index 6, 4 points
    10, // Precision 5: Starts at index 10, 5 points
    15, // Precision 6: Starts at index 15, 6 points
};

class GaussianQuadrature_1d {
public:

    virtual ~GaussianQuadrature_1d() {
    }

    /* Integrates function f(x) in the limits a to b, at the given precision */
    scalar integrate_function_1d(scalar a, scalar b, int precision);

    scalar integrate_function_1d_tri(scalar theta_max, scalar L_perp, scalar theta_star, int precision);

protected:
    virtual scalar f_1d(scalar r) = 0;
};
#endif
