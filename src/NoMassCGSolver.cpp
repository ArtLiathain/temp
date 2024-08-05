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

#include "NoMassCGSolver.h"

NoMassCGSolver::NoMassCGSolver() {
    num_rows = 0;
    num_nodes = 0;
    epsilon2 = 0;
    i_max = 0;
    preconditioner = NULL;
    r = NULL;
    p = NULL;
    z = NULL;
    q = NULL;
    f = NULL;
    V = NULL;
}

/* */
NoMassCGSolver::~NoMassCGSolver() {
    delete[] r;
    delete[] p;
    delete[] z;
    delete[] q;
    delete[] f;
    delete[] preconditioner;
    r = NULL;
    p = NULL;
    z = NULL;
    q = NULL;
    f = NULL;
    preconditioner = NULL;
    num_rows = 0;
    num_nodes = 0;
    epsilon2 = 0;
    i_max = 0;
    delete V;
    V = NULL;
}

/* */
int NoMassCGSolver::init(int num_nodes, int num_elements, mesh_node *node, tetra_element_linear *elem, SimulationParams *params, int num_pinned_nodes, int *pinned_nodes_list, set<int> bsite_pinned_node_list) {

    this->num_rows = 3 * num_nodes;
    this->num_nodes = num_nodes;
    this->epsilon2 = params->epsilon2;
    this->i_max = params->max_iterations_cg;
    this->one = 1;
    //printf("\t\t\tCalculating Sparsity Pattern for a 1st Order Viscosity Matrix\n");
    SparsityPattern sparsity_pattern_viscosity_matrix;
    sparsity_pattern_viscosity_matrix.init(num_rows);

    scalar *mem_loc;
    int n, ni, nj, ni_index, nj_index, ni_row, nj_row, i, j;
    matrix3 J;

    // Create a temporary lookup for checking if a node is 'pinned' or not.
    // if it is, then only a 1 on the diagonal corresponding to that node should
    // be placed (no off diagonal), effectively taking this node out of the equation
    // and therefore meaning the force on it should always be zero.
    vector<int> is_pinned(num_nodes);
    for (i = 0; i < num_nodes; i++) {
        is_pinned[i] = 0;
    }
    for (i = 0; i < num_pinned_nodes; i++) {
        is_pinned[pinned_nodes_list[i]] = 1;
    }
    for(set<int>::iterator it = bsite_pinned_node_list.begin(); it != bsite_pinned_node_list.end(); ++it) {
        is_pinned[*it] = 1;
    }

    for (n = 0; n < num_elements; n++) {
        elem[n].calculate_jacobian(J);
        elem[n].calc_shape_function_derivatives_and_volume(J);
        elem[n].create_viscosity_matrix();
        for (ni = 0; ni < 10; ++ni) {
            for (nj = 0; nj < 10; ++nj) {
                ni_index = elem[n].n[ni]->index;
                nj_index = elem[n].n[nj]->index;
		ni_row = ni_index * 3;
		nj_row = nj_index * 3;
                for (i = 0; i < 3; ++i) {
                    for (j = 0; j < 3; ++j) {
			if(is_pinned[ni_index] == 0 && is_pinned[nj_index] == 0) {
		            if (ni < 4 && nj < 4) {
		                mem_loc = &elem[n].viscosity_matrix[ni + 4 * i][nj + 4 * j];
		                sparsity_pattern_viscosity_matrix.register_contribution(ni_row + i, nj_row + j, mem_loc);
		            } else {
		                if (ni == nj && i == j) {
		                    if (sparsity_pattern_viscosity_matrix.check_for_contribution(ni_row + i, nj_row + j) == false) {
		                        mem_loc = &one;
		                        sparsity_pattern_viscosity_matrix.register_contribution(ni_row + i, nj_row + j, mem_loc);
		                    }
		                }
		            }
			} else {
			    if (ni == nj && i == j) {
		                if (sparsity_pattern_viscosity_matrix.check_for_contribution(ni_row + i, nj_row + j) == false) {
		                    mem_loc = &one;
		                    sparsity_pattern_viscosity_matrix.register_contribution(ni_row + i, nj_row + j, mem_loc);
		                }
		            }
			}
                    }
                }
            }
        }
    }

    if (params->calc_stokes == 1) {
        for (ni = 0; ni < num_nodes; ++ni) {
	    if(is_pinned[ni] == 0) {
                for (nj = 0; nj < 3; ++nj) {
                    sparsity_pattern_viscosity_matrix.register_contribution(3 * ni + nj, 3 * ni + nj, &node[ni].stokes_drag);
                }
	    }
        }
    }

    // Creating fixed pattern viscosity matrix, but not entering values! Just making the pattern for now
    //printf("\t\t\tBuilding Sparsity Pattern for Viscosity Matrix\n");
    V = sparsity_pattern_viscosity_matrix.create_sparse_matrix();
    //V->build();
    //delete V;
    //return FFEA_ERROR;
    // create a preconditioner for solving in less iterations
    // Create the jacobi preconditioner matrix (diagonal)
    preconditioner = new(std::nothrow) scalar[num_rows];
    if (preconditioner == NULL) FFEA_ERROR_MESSG("Failed to allocate 'preconditioner' in NoMassCGSolver\n");


    // create the work vectors necessary for use by the conjugate gradient solver in 'solve'
    r = new(std::nothrow) vector3[num_nodes];
    p = new(std::nothrow) vector3[num_nodes];
    z = new(std::nothrow) vector3[num_nodes];
    q = new(std::nothrow) vector3[num_nodes];
    f = new(std::nothrow) vector3[num_nodes];
    if (r == NULL || p == NULL || z == NULL || q == NULL || f == NULL) FFEA_ERROR_MESSG(" Failed to create the work vectors necessary for NoMassCGSolver\n"); 

    return FFEA_OK;
}

/*  */
int NoMassCGSolver::solve(vector3* x) {
    // Complete the sparse viscosity matrix
    V->build();
    V->calc_inverse_diagonal(preconditioner);
    //V->print_dense_to_file(x);
    //exit(0);
    int i = 0;
    scalar delta_new, delta_old, pTq, alpha;
    delta_new = conjugate_gradient_residual_assume_x_zero(x);

    for (i = 0; i < i_max; i++) {
        pTq = get_alpha_denominator();
        alpha = delta_new / pTq;

        // Update solution and residual
        vec3_add_to_scaled(x, p, alpha, num_nodes);
        vec3_add_to_scaled(r, q, -alpha, num_nodes);

        // Once convergence is achieved, return
        if (residual2() < epsilon2) {
	    //cout << residual2() << " " << epsilon2 << endl;
	    //exit(0);
            //std::cout << "NoMassCG_solver: Convergence reached on iteration " << i << "\n"; // DEBUGGO
            return FFEA_OK;
        }
	//cout << residual2() << " " << epsilon2 << endl;
        delta_old = delta_new;
        delta_new = parallel_apply_preconditioner();
        vec3_scale_and_add(p, z, (delta_new / delta_old), num_nodes);
    }
    // If desired convergence was not reached in the set number of iterations...
    FFEA_ERROR_MESSG("Conjugate gradient solver: Could not converge after %d iterations.\n\tEither epsilon or max_iterations_cg are set too low, or something went wrong with the simulation.\n", i_max);
}

/* */
void NoMassCGSolver::print_matrices(vector3* force) {
    V->print_dense_to_file(force);
}

/* */
scalar NoMassCGSolver::conjugate_gradient_residual_assume_x_zero(vector3 *b) {

    scalar delta_new = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) shared(b) reduction(+:delta_new)
#endif
    for (int i = 0; i < num_nodes; i++) {
        r[i].x = b[i].x;
        r[i].y = b[i].y;
        r[i].z = b[i].z;
        f[i].x = b[i].x;
        f[i].y = b[i].y;
        f[i].z = b[i].z;
        b[i].x = 0;
        b[i].y = 0;
        b[i].z = 0;
        z[i].x = preconditioner[(3 * i)] * r[i].x;
        z[i].y = preconditioner[(3 * i) + 1] * r[i].y;
        z[i].z = preconditioner[(3 * i) + 2] * r[i].z;
        p[i].x = z[i].x;
        p[i].y = z[i].y;
        p[i].z = z[i].z;
        delta_new += r[i].x * z[i].x + r[i].y * z[i].y + r[i].z * z[i].z;
    }
    return delta_new;
}

/* */
scalar NoMassCGSolver::residual2() {
    scalar r2 = 0, f2 = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) reduction(+:r2, f2)
#endif
    for (int i = 0; i < num_nodes; i++) {
        r2 += r[i].x * r[i].x + r[i].y * r[i].y + r[i].z * r[i].z;
        f2 += f[i].x * f[i].x + f[i].y * f[i].y + f[i].z * f[i].z;
    }
    if (f2 == 0.0) {
        return 0.0;
    } else {
        return r2 / f2;
    }
}

/* */
scalar NoMassCGSolver::modx(vector3 *x) {
    scalar r2 = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) shared(x) reduction(+:r2)
#endif
    for (int i = 0; i < num_nodes; i++) {
        r2 += x[i].x * x[i].x + x[i].y * x[i].y + x[i].z * x[i].z;
    }
    return r2;
}

scalar NoMassCGSolver::get_alpha_denominator() {
    // A * p
    V->apply(p, q);
    scalar pTq = 0;

    // p^T * A * p
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) reduction(+:pTq)
#endif
    for (int i = 0; i < num_nodes; ++i) {
        pTq += p[i].x * q[i].x + p[i].y * q[i].y + p[i].z * q[i].z;
    }

    return pTq;
}

/* */
scalar NoMassCGSolver::parallel_apply_preconditioner() {
    scalar delta_new = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) reduction(+:delta_new)
#endif
    for (int i = 0; i < num_nodes; i++) {
        z[i].x = preconditioner[(3 * i)] * r[i].x;
        z[i].y = preconditioner[(3 * i) + 1] * r[i].y;
        z[i].z = preconditioner[(3 * i) + 2] * r[i].z;
        delta_new += r[i].x * z[i].x + r[i].y * z[i].y + r[i].z * z[i].z;
    }

    return delta_new;
}

/* */
void NoMassCGSolver::check(vector3 *x) {
    FILE *fout2;
    fout2 = fopen("/localhome/py09bh/output/nomass/cube_viscosity_no_mass.csv", "a");
    int i;
    double temp = 0, temp2 = 0;
    vector3 *temp_vec = new vector3[num_nodes];
    V->apply(x, temp_vec);
    for (i = 0; i < num_nodes; ++i) {
        temp += x[i].x * temp_vec[i].x + x[i].y * temp_vec[i].y + x[i].z * temp_vec[i].z;
        temp2 += x[i].x * f[i].x + x[i].y * f[i].y + x[i].z * f[i].z;
    }
    fprintf(fout2, "%e,%e\n", temp2, fabs(temp - temp2));
    fclose(fout2);

}

/* */



