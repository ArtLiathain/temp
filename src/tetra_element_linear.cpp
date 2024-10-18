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
 * tetra_element_linear.cpp
 */

#include "tetra_element_linear.h"

tetra_element_linear::tetra_element_linear() {
    rho = 0;
    A = 0;
    B = 0;
    G = 0;
    E = 0;
    dielectric = 0;
    // Not sure equivalent now switched to std::array
    //for (int i = 0; i < NUM_NODES_QUADRATIC_TET; i++) {
    //    n[i] = nullptr;
    //}

    //			node_phi[0] = 0; node_phi[1] = 0; node_phi[2] = 0; node_phi[3] = 0;
    vol_0 = 0;
    vol = 0;
    mat3_set_identity(F_ij);
    internal_stress_mag = 0;
    mat3_set_zero(J_inv_0);
    mat12_set_zero(viscosity_matrix);
    zero_force();
    last_det = 0;
    daddy_blob = nullptr;
}

/*
 * Get the memory location of the specified element of K_alpha
 */
scalar * tetra_element_linear::get_K_alpha_element_mem_loc(int ni, int nj) {
    return K_alpha.get_K_alpha_mem_loc(ni, nj);
}

/* Calc the diffusion matrix for this element */
void tetra_element_linear::calculate_K_alpha() {
    // Build the poisson diffusion matrix corresponding to this element
    K_alpha.build(n, dielectric);
}

void tetra_element_linear::construct_element_mass_matrix(MassMatrixQuadratic &M_alpha) {
    // Build the element mass matrix corresponding to this element
    M_alpha.build(n);
}

void tetra_element_linear::construct_element_mass_matrix(MassMatrixLinear &M_alpha) {
    // Build the element mass matrix corresponding to this element
    M_alpha.build(rho, vol_0);
}

void tetra_element_linear::add_K_alpha(scalar *K, int num_nodes) {
    for (int i = 0; i < NUM_NODES_QUADRATIC_TET; i++) {
        for (int j = 0; j < NUM_NODES_QUADRATIC_TET; j++) {
            K[n[i]->index * num_nodes + n[j]->index] += K_alpha.get_K_alpha_value(i, j);
        }
    }
}

/* Returns the gradient of the potential at the given (s,t,u) position in the element */
void tetra_element_linear::get_grad_phi_at_stu(arr3 &grad_phi, scalar s, scalar t, scalar u) {
    std::array<arr3, NUM_NODES_QUADRATIC_TET> grad_psi = {};

    SecondOrderFunctions::abcd3x3 J_coeff = {};
    SecondOrderFunctions::calc_jacobian_column_coefficients(n, J_coeff);
    vector9 J_inv;
    scalar det_J = SecondOrderFunctions::calc_det_J(J_coeff, s, t, u, J_inv);
    SecondOrderFunctions::calc_grad_psi(grad_psi, s, t, u, J_inv);

    grad_phi[0] = 0;
    grad_phi[1] = 0;
    grad_phi[2] = 0;
    for (int i = 0; i < NUM_NODES_QUADRATIC_TET; i++) {
        grad_phi[0] += grad_psi[i][0] * n[i]->phi;
        grad_phi[1] += grad_psi[i][1] * n[i]->phi;
        grad_phi[2] += grad_psi[i][2] * n[i]->phi;
    }

    grad_phi[0] *= det_J;
    grad_phi[1] *= det_J;
    grad_phi[2] *= det_J;
}

/* Calculates the force on each node of the element due to the electrostatic potential gradient there */
void tetra_element_linear::calculate_electrostatic_forces() {
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
    std::array<arr3, NUM_NODES_QUADRATIC_TET> grad_psi = {};
    std::array<scalar, NUM_NODES_QUADRATIC_TET> psi = {};
    std::array<arr3, NUM_NODES_QUADRATIC_TET> force = {};
    arr3 grad_phi_here;
    

    SecondOrderFunctions::calc_jacobian_column_coefficients(n, J_coeff);

    vector9 J_inv = {};

    for (int i = 0; i < NUM_TET_GAUSS_QUAD_POINTS; i++) {
        scalar det_J = SecondOrderFunctions::calc_det_J(J_coeff, gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], J_inv);
        SecondOrderFunctions::calc_grad_psi(grad_psi, gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], J_inv);
        SecondOrderFunctions::calc_psi(psi, gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2]);

        grad_phi_here[0] = 0;
        grad_phi_here[1] = 0;
        grad_phi_here[2] = 0;
        for (int j = 0; j < NUM_NODES_QUADRATIC_TET; j++) {
            grad_phi_here[0] += grad_psi[j][0] * n[j]->phi * det_J;
            grad_phi_here[1] += grad_psi[j][1] * n[j]->phi * det_J;
            grad_phi_here[2] += grad_psi[j][2] * n[j]->phi * det_J;
            //					printf("%d >> %e (%e %e %e)\n", i, det_J, grad_psi[j][0], grad_psi[j][1], grad_psi[j][2]);
        }

        //				printf("%d :: (%e, %e, %e)\n", i, grad_phi_here[0], grad_phi_here[1], grad_phi_here[2]);

        for (int j = 0; j < NUM_NODES_QUADRATIC_TET; j++) {
            //					force[j][0] -= gauss_points[i].W * grad_phi_here[0] * psi[j] * n[j]->rho;
            //					force[j][1] -= gauss_points[i].W * grad_phi_here[1] * psi[j] * n[j]->rho;
            //					force[j][2] -= gauss_points[i].W * grad_phi_here[2] * psi[j] * n[j]->rho;
            //					printf("%d == %e\n", i, n[j]->rho);
            force[j][0] -= gauss_points[i].W * grad_phi_here[0] * n[j]->rho;
            force[j][1] -= gauss_points[i].W * grad_phi_here[1] * n[j]->rho;
            force[j][2] -= gauss_points[i].W * grad_phi_here[2] * n[j]->rho;
        }
    }

    for (int i = 0; i < NUM_NODES_QUADRATIC_TET; i++) {
        add_force_to_node(i, force[i]);
    }

}

/* Calculates the Jacobian matrix for this element */
void tetra_element_linear::calculate_jacobian(matrix3 &J) {
    J[0][0] = n[1]->pos[0] - n[0]->pos[0];
    J[0][1] = n[1]->pos[1] - n[0]->pos[1];
    J[0][2] = n[1]->pos[2] - n[0]->pos[2];

    J[1][0] = n[2]->pos[0] - n[0]->pos[0];
    J[1][1] = n[2]->pos[1] - n[0]->pos[1];
    J[1][2] = n[2]->pos[2] - n[0]->pos[2];

    J[2][0] = n[3]->pos[0] - n[0]->pos[0];
    J[2][1] = n[3]->pos[1] - n[0]->pos[1];
    J[2][2] = n[3]->pos[2] - n[0]->pos[2];
}

/*
 * Inverts the given jacobian matrix J, using this to calculate the derivatives
 * of the shape functions which are stored in dpsi. This is an array
 * of all 12 derivatives, in the following order:
 *		dpsi    =	[ d(psi)_1/dx ]
 *				[ d(psi)_2/dx ]
 *				[ d(psi)_3/dx ]
 *				[ d(psi)_4/dx ]
 *				[ d(psi)_1/dy ]
 *				[ d(psi)_2/dy ]
 *				[     ...     ]
 *				[ d(psi)_4/dz ]
 *
 * Function also (as a by-product of the inversion) calculates the volume of the
 * element whose jacobian this is, which is stored in 'vol'.
 * @return True if the element has inverted
 */
bool tetra_element_linear::calc_shape_function_derivatives_and_volume(matrix3 &J) {
    scalar det;

    // Calculate shape function derivs from inverse jacobian directly into dpsi[]

    /*
                            dpsi[1] = J[2][2]*J[1][1] - J[2][1]*J[1][2];
                            dpsi[2] = J[2][1]*J[0][2] - J[2][2]*J[0][1];
                            dpsi[3] = J[1][2]*J[0][1] - J[1][1]*J[0][2];

                            dpsi[5] = J[2][0]*J[1][2] - J[2][2]*J[1][0];
                            dpsi[6] = J[2][2]*J[0][0] - J[2][0]*J[0][2];
                            dpsi[7] = J[1][0]*J[0][2] - J[1][2]*J[0][0];

                            dpsi[9] = J[2][1]*J[1][0] - J[2][0]*J[1][1];
                            dpsi[10] = J[2][0]*J[0][1] - J[2][1]*J[0][0];
                            dpsi[11] = J[1][1]*J[0][0] - J[1][0]*J[0][1];
     */
    dpsi[DPSI2_DX] = J[2][2] * J[1][1] - J[2][1] * J[1][2];
    dpsi[DPSI3_DX] = J[2][1] * J[0][2] - J[2][2] * J[0][1];
    dpsi[DPSI4_DX] = J[1][2] * J[0][1] - J[1][1] * J[0][2];

    dpsi[DPSI2_DY] = J[2][0] * J[1][2] - J[2][2] * J[1][0];
    dpsi[DPSI3_DY] = J[2][2] * J[0][0] - J[2][0] * J[0][2];
    dpsi[DPSI4_DY] = J[1][0] * J[0][2] - J[1][2] * J[0][0];

    dpsi[DPSI2_DZ] = J[2][1] * J[1][0] - J[2][0] * J[1][1];
    dpsi[DPSI3_DZ] = J[2][0] * J[0][1] - J[2][1] * J[0][0];
    dpsi[DPSI4_DZ] = J[1][1] * J[0][0] - J[1][0] * J[0][1];

    // Determinant
    det = J[0][0] * dpsi[DPSI2_DX] + J[1][0] * dpsi[DPSI3_DX] + J[2][0] * dpsi[DPSI4_DX];

    // Check if element has inverted itself (determinant changed sign)
 /*   if(index == 681) {
	printf("J Dets: Last = %e Now = %e\n", last_det, det);
    } */
    if (last_det * det < 0) {
        return true;
    }
    last_det = det;

    // Calculate volume of element from the determinant
    // and the volume of the "parent" element in local coords
    // e.g 1/6 for a right angled tetrahedron
    vol = (1.0 / 6.0) * fabs(det);

    // Divide by determinant to complete matrix inversion
    det = 1.0 / det;
    dpsi[1] *= det;
    dpsi[2] *= det;
    dpsi[3] *= det;
    dpsi[5] *= det;
    dpsi[6] *= det;
    dpsi[7] *= det;
    dpsi[9] *= det;
    dpsi[10] *= det;
    dpsi[11] *= det;

    // calculate the derivatives of shape function psi_1 w.r.t. x, y and z
    /*
                            dpsi[0] = -(dpsi[1] + dpsi[2] + dpsi[3]);
                            dpsi[4] = -(dpsi[5] + dpsi[6] + dpsi[7]);
                            dpsi[8] = -(dpsi[9] + dpsi[10] + dpsi[11]);
     */
    dpsi[DPSI1_DX] = -(dpsi[DPSI2_DX] + dpsi[DPSI3_DX] + dpsi[DPSI4_DX]);
    dpsi[DPSI1_DY] = -(dpsi[DPSI2_DY] + dpsi[DPSI3_DY] + dpsi[DPSI4_DY]);
    dpsi[DPSI1_DZ] = -(dpsi[DPSI2_DZ] + dpsi[DPSI3_DZ] + dpsi[DPSI4_DZ]);

    return false;
}

/*
 * Builds the viscosity matrix from the shape function derivatives, the shear and bulk
 * viscosity constants, and the element volume
 */
void tetra_element_linear::create_viscosity_matrix() {
    int i, j;
    matrix4 K;

    // Construct submatrices on the diagonal
    BULK_VISCOUS_SUBMATRIX_DIAG(viscosity_matrix, 0)
    BULK_VISCOUS_SUBMATRIX_DIAG(viscosity_matrix, 4)
    BULK_VISCOUS_SUBMATRIX_DIAG(viscosity_matrix, 8)

    // Construct submatrices off diagonal
    BULK_VISCOUS_SUBMATRIX_OFFDIAG(viscosity_matrix, 0, 4)
    BULK_VISCOUS_SUBMATRIX_OFFDIAG(viscosity_matrix, 0, 8)
    BULK_VISCOUS_SUBMATRIX_OFFDIAG(viscosity_matrix, 4, 8)

    // Create the diffusion matrix and add it to the viscosity matrix in 3 upper
    // triangular blocks along the diagonal
    calc_del2_matrix();
    add_diffusion_matrix(viscosity_matrix);

    // Multiply all these values (currently only in upper half of matrix) by the element volume
    for (i = 0; i < 12; i++)
        for (j = 0; j <= i; j++)
            viscosity_matrix[j][i] *= vol;

    // Viscosity matrix is symmetric, so no need to recalculate entries.
    // Simply 'mirror image' the matrix back into itself
    for (i = 1; i < 12; i++)
        for (j = 0; j < i; j++)
            viscosity_matrix[i][j] = viscosity_matrix[j][i];
}

/*
 *
 */
void tetra_element_linear::print_structural_details() {

	printf("Element %d\n\n", index);
	printf("\tVolume = %5.2f angstroms\n", calc_volume() * 4.913);	// Making it into angstroms
	printf("\n");
}

/*
 *
 */
scalar tetra_element_linear::calc_volume() {

	matrix3 J = {};
	calculate_jacobian(J);
        calc_shape_function_derivatives_and_volume(J);
	return vol;
}

/*
 *
 */
void tetra_element_linear::calc_elastic_force_vector(vector12 &F) {
	
	matrix3 J = {}, stress = {};
	mat3_set_zero(stress);
	vec12_set_zero(F);

	calculate_jacobian(J);
        calc_shape_function_derivatives_and_volume(J);
	add_shear_elastic_stress(J, stress);
	add_bulk_elastic_stress(stress);
	apply_stress_tensor(stress, F);
}

/*
 *
 */
void tetra_element_linear::calc_deformation(matrix3 &J) {
	// Reset gradient deformation to zero
   	mat3_set_zero(F_ij);

	// F_ij transpose is given by the current jacobian times the rest state jacobian inverse.
    	mat3_mult_both_transposed(J, J_inv_0, F_ij);
}

void tetra_element_linear::add_shear_elastic_stress(matrix3 &J, matrix3 &stress) {

    // Reset gradient deformation to zero
    mat3_set_zero(F_ij);

    // F_ij transpose is given by the current jacobian times the rest state jacobian inverse.
    mat3_mult_both_transposed(J, J_inv_0, F_ij);

    // Now multiply F_ij with its own transpose to get F_ik_kjT
    mat3_mult_transpose(F_ij, F_ij, stress);

    // Scale matrix by the appropriate factors
    mat3_scale(stress, (G * vol_0 / vol));

    // Subtract Identity to make deformation tensor (almost) traceless
    stress[0][0] -= G;
    stress[1][1] -= G;
    stress[2][2] -= G;
}

/*
 *
 */
void tetra_element_linear::add_bulk_elastic_stress(matrix3 &stress) {

    scalar c_2 = E - G * 2.0 / 3.0;
    scalar c = G * (1.0 - (vol_0 / vol)) + 0.5 * c_2 * ((vol / vol_0) - (vol_0 / vol));
    stress[0][0] += c;
    stress[1][1] += c;
    stress[2][2] += c;
}

/*
 * Given the shape function derivatives, the element volume and a random number generator, this
 * function calculates the fluctuating stress tensor, generating a stochastic change in the
 * nodal velocities for the element under consideration. This function will add its contribution
 * to the given 12-vector du.
 *
 */
void tetra_element_linear::add_fluctuating_stress(const SimulationParams &params, std::shared_ptr<std::vector<RngStream>> &rng, matrix3 stress, int thread_id) {
    scalar c = sqrt((24 * params.kT) / (vol * params.dt));

    // Bulk fluctuation term
    scalar bf = sqrt(B) * RAND(-.5, .5);

    // Diagonal terms
    stress[0][0] += c * (sqrt(2 * A) * RAND(-.5, .5) + bf);
    stress[1][1] += c * (sqrt(2 * A) * RAND(-.5, .5) + bf);
    stress[2][2] += c * (sqrt(2 * A) * RAND(-.5, .5) + bf);

    // Off diagonal terms (note that stress matrix is symmetric)
    stress[0][1] += c * sqrt(A) * RAND(-.5, .5);
    stress[0][2] += c * sqrt(A) * RAND(-.5, .5);
    stress[1][2] += c * sqrt(A) * RAND(-.5, .5);

    stress[1][0] = stress[0][1];
    stress[2][0] = stress[0][2];
    stress[2][1] = stress[1][2];
}

/*
 * Applies the given stress tensor to the shape function derivatives to get the contribution to du
 */
void tetra_element_linear::apply_stress_tensor(matrix3 &stress, vector12 &du) {
    int i;
    for (i = 0; i < 3; i++) {
        du[4 * i] += vol * (dpsi[0] * stress[i][0] + dpsi[4] * stress[i][1] + dpsi[8] * stress[i][2]);
        du[4 * i + 1] += vol * (dpsi[1] * stress[i][0] + dpsi[5] * stress[i][1] + dpsi[9] * stress[i][2]);
        du[4 * i + 2] += vol * (dpsi[2] * stress[i][0] + dpsi[6] * stress[i][1] + dpsi[10] * stress[i][2]);
        du[4 * i + 3] += vol * (dpsi[3] * stress[i][0] + dpsi[7] * stress[i][1] + dpsi[11] * stress[i][2]);
    }
}

/*
 * Sets the given 12-vector to the velocities of this element's four nodes,
 */
void tetra_element_linear::get_element_velocity_vector(vector12 &v) {
    v[0] = n[0]->vel[0];
    v[1] = n[1]->vel[0];
    v[2] = n[2]->vel[0];
    v[3] = n[3]->vel[0];

    v[4] = n[0]->vel[1];
    v[5] = n[1]->vel[1];
    v[6] = n[2]->vel[1];
    v[7] = n[3]->vel[1];

    v[8] = n[0]->vel[2];
    v[9] = n[1]->vel[2];
    v[10] = n[2]->vel[2];
    v[11] = n[3]->vel[2];
}

/*
 * Add this element's nodal forces to those given in the force 12-vector
 */
void tetra_element_linear::add_element_force_vector(vector12 &force) {
    node_force[0][0] -= force[0];
    node_force[1][0] -= force[1];
    node_force[2][0] -= force[2];
    node_force[3][0] -= force[3];

    node_force[0][1] -= force[4];
    node_force[1][1] -= force[5];
    node_force[2][1] -= force[6];
    node_force[3][1] -= force[7];

    node_force[0][2] -= force[8];
    node_force[1][2] -= force[9];
    node_force[2][2] -= force[10];
    node_force[3][2] -= force[11];
}

/* Add given force to the specified node of this element */
void tetra_element_linear::add_force_to_node(int i, arr3 &f) {
    node_force[i][0] += f[0];
    node_force[i][1] += f[1];
    node_force[i][2] += f[2];
}

/* A roundabout and inefficient way of working out what node (from 0 to 9) this index corresponds to */
int tetra_element_linear::what_node_is_this(int index) {
    for (int i = 0; i < NUM_NODES_QUADRATIC_TET; i++) {
        if (n[i]->index == index) {
            return i;
        }
    }
    throw FFEAException("Specified node index does not belong to this element.");
}

void tetra_element_linear::print() {
    int i;
    for (i = 0; i < NUM_NODES_QUADRATIC_TET; i++) {
        printf("Node %d:\n", i);
        n[i]->print();
        printf("node_force: %e %e %e\n", node_force[i][0], node_force[i][1], node_force[i][2]);
        printf("volume: %e\n", vol);
    }
}

void tetra_element_linear::print_viscosity_matrix() {
	for(int i = 0; i < 12; ++i) {
		for(int j = 0; j < 12; ++j) {
			cout << viscosity_matrix[i][j] << " ";
		}
		cout << endl;
	}
}

/*
 * Applies the mass matrix (for a linear tetrahedral element of density rho and equilibrium volume vol_0)
 * to the force vector du to get the correct distribution of force between nodes.
 */
void tetra_element_linear::apply_element_mass_matrix(vector12 &du) {
    int i, j, k;

    // mass matrix
    const matrix4 M = {
        std::array{.1, .05, .05, .05},
        std::array{.05, .1, .05, .05},
        std::array{.05, .05, .1, .05},
        std::array{.05, .05, .05, .1}};

    vector12 temp_du = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    // Apply mass matrix to each of the three subvectors
    for (k = 0; k < 12; k += 4)
        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
                temp_du[i + k] += M[i][j] * du[j + k];

    for (i = 0; i < 12; i++) du[i] = temp_du[i] * (vol_0 * rho);
}

void tetra_element_linear::volume_coord_to_xyz(scalar eta0, scalar eta1, scalar eta2, scalar eta3, arr3 &r) {
    for(int i = 0; i < 3; ++i)
        r[i] = eta0 * n[0]->pos[i] + eta1 * n[1]->pos[i] + eta2 * n[2]->pos[i] + eta3 * n[3]->pos[i];
}

void tetra_element_linear::zero_force() {
    for (int i = 0; i < NUM_NODES_QUADRATIC_TET; i++) {
        arr3_set_zero(node_force[i]);
    }
}

void tetra_element_linear::linearise_element() {
    n[4]->pos[0] = .5 * (n[0]->pos[0] + n[1]->pos[0]);
    n[4]->pos[1] = .5 * (n[0]->pos[1] + n[1]->pos[1]);
    n[4]->pos[2] = .5 * (n[0]->pos[2] + n[1]->pos[2]);

    n[5]->pos[0] = .5 * (n[0]->pos[0] + n[2]->pos[0]);
    n[5]->pos[1] = .5 * (n[0]->pos[1] + n[2]->pos[1]);
    n[5]->pos[2] = .5 * (n[0]->pos[2] + n[2]->pos[2]);

    n[6]->pos[0] = .5 * (n[0]->pos[0] + n[3]->pos[0]);
    n[6]->pos[1] = .5 * (n[0]->pos[1] + n[3]->pos[1]);
    n[6]->pos[2] = .5 * (n[0]->pos[2] + n[3]->pos[2]);

    n[7]->pos[0] = .5 * (n[1]->pos[0] + n[2]->pos[0]);
    n[7]->pos[1] = .5 * (n[1]->pos[1] + n[2]->pos[1]);
    n[7]->pos[2] = .5 * (n[1]->pos[2] + n[2]->pos[2]);

    n[8]->pos[0] = .5 * (n[1]->pos[0] + n[3]->pos[0]);
    n[8]->pos[1] = .5 * (n[1]->pos[1] + n[3]->pos[1]);
    n[8]->pos[2] = .5 * (n[1]->pos[2] + n[3]->pos[2]);

    n[9]->pos[0] = .5 * (n[2]->pos[0] + n[3]->pos[0]);
    n[9]->pos[1] = .5 * (n[2]->pos[1] + n[3]->pos[1]);
    n[9]->pos[2] = .5 * (n[2]->pos[2] + n[3]->pos[2]);
}

void tetra_element_linear::calc_centroid() {
    centroid[0] = .25 * (n[0]->pos[0] + n[1]->pos[0] + n[2]->pos[0] + n[3]->pos[0]);
    centroid[1] = .25 * (n[0]->pos[1] + n[1]->pos[1] + n[2]->pos[1] + n[3]->pos[1]);
    centroid[2] = .25 * (n[0]->pos[2] + n[1]->pos[2] + n[2]->pos[2] + n[3]->pos[2]);
}

void tetra_element_linear::calc_del2_matrix() {
    del2.u00 = (dpsi[0] * dpsi[0] + dpsi[4] * dpsi[4] + dpsi[8] * dpsi[8]);
    del2.u01 = (dpsi[0] * dpsi[1] + dpsi[4] * dpsi[5] + dpsi[8] * dpsi[9]);
    del2.u02 = (dpsi[0] * dpsi[2] + dpsi[4] * dpsi[6] + dpsi[8] * dpsi[10]);
    del2.u03 = (dpsi[0] * dpsi[3] + dpsi[4] * dpsi[7] + dpsi[8] * dpsi[11]);

    del2.u11 = (dpsi[1] * dpsi[1] + dpsi[5] * dpsi[5] + dpsi[9] * dpsi[9]);
    del2.u12 = (dpsi[1] * dpsi[2] + dpsi[5] * dpsi[6] + dpsi[9] * dpsi[10]);
    del2.u13 = (dpsi[1] * dpsi[3] + dpsi[5] * dpsi[7] + dpsi[9] * dpsi[11]);

    del2.u22 = (dpsi[2] * dpsi[2] + dpsi[6] * dpsi[6] + dpsi[10] * dpsi[10]);
    del2.u23 = (dpsi[2] * dpsi[3] + dpsi[6] * dpsi[7] + dpsi[10] * dpsi[11]);

    del2.u33 = (dpsi[3] * dpsi[3] + dpsi[7] * dpsi[7] + dpsi[11] * dpsi[11]);

    //			printf("del2:\n");
    //			printf("%e %e %e %e\n", del2.u00, del2.u01, del2.u02, del2.u03);
    //			printf("%e %e %e %e\n", del2.u01, del2.u11, del2.u12, del2.u13);
    //			printf("%e %e %e %e\n", del2.u02, del2.u12, del2.u22, del2.u23);
    //			printf("%e %e %e %e\n", del2.u03, del2.u13, del2.u23, del2.u33);
}

void tetra_element_linear::add_diffusion_matrix(matrix12 &V) {
    int i;

    // Drop each upper triangle of this diffusion matrix along the block diagonal
    // of the viscosity matrix
    for (i = 0; i < 3; i++) {
        V[i * 4 + 0][i * 4 + 0] += del2.u00 * A;
        V[i * 4 + 0][i * 4 + 1] += del2.u01 * A;
        V[i * 4 + 0][i * 4 + 2] += del2.u02 * A;
        V[i * 4 + 0][i * 4 + 3] += del2.u03 * A;
        V[i * 4 + 1][i * 4 + 1] += del2.u11 * A;
        V[i * 4 + 1][i * 4 + 2] += del2.u12 * A;
        V[i * 4 + 1][i * 4 + 3] += del2.u13 * A;
        V[i * 4 + 2][i * 4 + 2] += del2.u22 * A;
        V[i * 4 + 2][i * 4 + 3] += del2.u23 * A;
        V[i * 4 + 3][i * 4 + 3] += del2.u33 * A;
    }
}

/** Returns the opposite node for a 2nd order face,
 *    e. g., returns 3 for face [0,4,5].
 *  The 2nd order nodes were created at setup time, in:
 *    FFEA_initialise/FFEA_convert_from_volume/convert_tetrahedra_linear_to_quadratic.py
 */
int tetra_element_linear::get_opposite_node(int n1, int n2, int n3) {
  
   // get a hash for the face:
   int hash = n1 * n1 + n2 * n2 + n3 * n3;

   if (hash == 114 || hash == 134 || hash == 154 || hash == 194) {
      return 0;
   } else
   if (hash == 110 || hash == 142 || hash == 126 || hash == 61) {
      return 1;
   } else
   if (hash == 116 || hash == 81 || hash == 109 || hash == 52) {
      return 2;
   } else
   if (hash == 66 || hash == 41 || hash == 78 || hash == 90) {
      return 3;
   } else return -1;

   
   /* 
   map<int,int> opposite;
   opposite[41]  = 3; // [0, 4, 5]
   opposite[52]  = 2; // [52, [0, 4, 6]]   # 2
   opposite[61]  = 1; // [61, [0, 5, 6]]   # 1
   opposite[66]  = 3; // [66, [1, 4, 7]]   # 3
   opposite[81]  = 2; // [81, [1, 4, 8]]   # 2
   opposite[114] = 0; // [114, [1, 7, 8]]  # 0
   opposite[110] = 1; // [110, [2, 5, 9]]  # 1
   opposite[78]  = 3; // [78, [2, 5, 7]]   # 3
   opposite[134] = 0; // [134, [2, 7, 9]]  # 0
   opposite[109] = 2; // [109, [3, 6, 8]]  # 2
   opposite[126] = 1; // [126, [3, 6, 9]]  # 1
   opposite[154] = 0; // [154, [3, 8, 9]]  # 0
   opposite[116] = 2; // [116, [4, 6, 8]]  # 2
   opposite[90]  = 3; // [90, [4, 5, 7]]   # 3
   opposite[142] = 1; // [142, [5, 6, 9]]  # 1
   opposite[194] = 0; // [194, [7, 8, 9]]  # 0

   // check whether this hash exists and:
   //   either fail:
   if (opposite.count(hash) != 1) return -1;
   //   or return the opposite node
   else return opposite[hash];
   */ 

}

scalar tetra_element_linear::length_of_longest_edge() {

   scalar d2 = 0;
   scalar di2 = 0;
   for (int i=0; i<NUM_NODES_LINEAR_TET; i++) { // loop over the linear nodes.
      for (int j=i+1; j<NUM_NODES_LINEAR_TET; j++) {
         di2 = (n[i]->pos[0] - n[j]->pos[0])*(n[i]->pos[0] - n[j]->pos[0]); 
         di2 += (n[i]->pos[1] - n[j]->pos[1])*(n[i]->pos[1] - n[j]->pos[1]);
         di2 += (n[i]->pos[2] - n[j]->pos[2])*(n[i]->pos[2] - n[j]->pos[2]);
         if (di2 > d2) d2 = di2;
      }
   }

   return sqrt(d2);

}

