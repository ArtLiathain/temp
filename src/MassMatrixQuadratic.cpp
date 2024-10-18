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

#include "MassMatrixQuadratic.h"

MassMatrixQuadratic::MassMatrixQuadratic() {
    zero();
}

scalar * MassMatrixQuadratic::get_M_alpha_mem_loc(int i, int j) {
    // Poisson matrix is symmetric, so convert any request for an upper triangular element into its
    // corresponding (equivalent) lower triangular element
    if (i < j) {
        int temp = i;
        i = j;
        j = temp;
    }

    // Get the index of the element corresponding to the position i,j in our lower triangular K_alpha matrix
    int c = (i * (i + 1)) / 2 + j;

    if (c < 0 || c > 54) {
        return nullptr;
    }

    // Return a pointer to this memory location
    return &(M_alpha[c]);
}

void MassMatrixQuadratic::build(std::array<mesh_node*, NUM_NODES_QUADRATIC_TET> &n) {
    struct tetrahedron_gauss_point gauss_points[NUM_TET_GAUSS_QUAD_POINTS] ={
        // Weight, eta1, eta2, eta3, eta4
        {0.317460317460317450e-2,
            {0.5, 0.5, 0.0, 0.0}},
        {0.317460317460317450e-2,
            {0.5, 0.0, 0.5, 0.0}},
        {0.317460317460317450e-2,
            {0.5, 0.0, 0.0, 0.5}},
        {0.317460317460317450e-2,
            {0.0, 0.5, 0.5, 0.0}},
        {0.317460317460317450e-2,
            {0.0, 0.0, 0.5, 0.5}},
        {0.317460317460317450e-2,
            {0.0, 0.5, 0.0, 0.5}},
        {0.147649707904967828e-1,
            {0.100526765225204467, 0.100526765225204467, 0.100526765225204467, 0.698419704324386603}},
        {0.147649707904967828e-1,
            {0.100526765225204467, 0.100526765225204467, 0.698419704324386603, 0.100526765225204467}},
        {0.147649707904967828e-1,
            {0.100526765225204467, 0.698419704324386603, 0.100526765225204467, 0.100526765225204467}},
        {0.147649707904967828e-1,
            {0.698419704324386603, 0.100526765225204467, 0.100526765225204467, 0.100526765225204467}},
        {0.221397911142651221e-1,
            {0.314372873493192195, 0.314372873493192195, 0.314372873493192195, 0.568813795204234229e-1}},
        {0.221397911142651221e-1,
            {0.314372873493192195, 0.314372873493192195, 0.568813795204234229e-1, 0.314372873493192195}},
        {0.221397911142651221e-1,
            {0.314372873493192195, 0.568813795204234229e-1, 0.314372873493192195, 0.314372873493192195}},
        {0.221397911142651221e-1,
            {0.568813795204234229e-1, 0.314372873493192195, 0.314372873493192195, 0.314372873493192195}}
    };

    SecondOrderFunctions::abcd3x3 J_coeff = {};
    std::array<scalar, NUM_NODES_QUADRATIC_TET> psi = {};

    SecondOrderFunctions::calc_jacobian_column_coefficients(n, J_coeff);

    zero();

    vector9 J_inv;

    for (int i = 0; i < NUM_TET_GAUSS_QUAD_POINTS; i++) {
        scalar det_J = SecondOrderFunctions::calc_det_J(J_coeff, gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], J_inv);
        SecondOrderFunctions::calc_psi(psi, gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2]);
        //				if(i == 0) printf("MassMatrix: det_J = %e\n", det_J);
        add_psi_dot_products(psi, det_J, gauss_points[i].W);
    }


    /*
                            printf("MASSMATRIXQUADRATIC:\n");
                            scalar M[10][10];
                            int c = 0;
                            for(int j = 0; j < 10; j++) {
                                    for(int i = 0; i <= j; i++) {
                                            M[i][j] = M_alpha[c];
                                            M[j][i] = M_alpha[c];
                                            c++;
                                    }
                            }

                            for(int i = 0; i < 10; i++) {
                                    for(int j = 0; j < 10; j++) {
                                            printf("%+e ", M[i][j]);
                                    }
                                    printf("\n");
                            }
                            printf("\n");
     */
}

scalar MassMatrixQuadratic::get_M_alpha_value(int i, int j) {
    // Poisson matrix is symmetric, so convert any request for an upper triangular element into its
    // corresponding (equivalent) lower triangular element
    if (i < j) {
        int temp = i;
        i = j;
        j = temp;
    }

    // Get the index of the element corresponding to the position i,j in our lower triangular K_alpha matrix
    int c = (i * (i + 1)) / 2 + j;

    if (c < 0 || c > 54) {
        return -897000;
    }

    // Return a pointer to this memory location
    return M_alpha[c];
}

void MassMatrixQuadratic::add_psi_dot_products(std::array<scalar, NUM_NODES_QUADRATIC_TET> &psi, scalar det_J, scalar weight) {
    int c = 0;
    for (int i = 0; i < NUM_NODES_QUADRATIC_TET; i++) {
        for (int j = 0; j <= i; j++) {
            M_alpha[c] += det_J * weight * psi[i] * psi[j];
            c++;
        }
    }
}

void MassMatrixQuadratic::zero() {
    for (int i = 0; i < NUM_ELEMENTS_LOWER_TRIANGULAR_10X10; i++) {
        M_alpha[i] = 0;
    }
}

void MassMatrixQuadratic::print_details() {
    // Printing total mass
    int c = 0;
    scalar tot = 0.0;
    for(int i = 0; i < 10; ++i) {
	for(int j = i; j < 10; ++j) {
	    if(i != j) {
		tot += M_alpha[c];
	    } else {
		tot += M_alpha[c];
	    }
	    c++;
	}
    }
    printf("Total Mass = %e\n", tot * mesoDimensions::mass);
}
