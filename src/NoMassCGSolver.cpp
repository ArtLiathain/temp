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
    preconditioner = {};
    r = {};
    p = {};
    z = {};
    q = {};
    f = {};
    V = nullptr;
}

/* */
NoMassCGSolver::~NoMassCGSolver() {
    r.clear();
    p.clear();
    z.clear();
    q.clear();
    f.clear();
    preconditioner.clear();
    num_rows = 0;
    num_nodes = 0;
    epsilon2 = 0;
    i_max = 0;
    V.reset();
}

/* */
void NoMassCGSolver::init(std::vector<mesh_node> &node, std::vector<tetra_element_linear> &elem, const SimulationParams &params, const std::vector<int> &pinned_nodes_list, const set<int> &bsite_pinned_node_list) {
    this->num_rows = 3 * node.size();
    this->num_nodes = node.size();
    this->epsilon2 = params.epsilon2;
    this->i_max = params.max_iterations_cg;
    this->one = 1;
    //printf("\t\t\tCalculating Sparsity Pattern for a 1st Order Viscosity Matrix\n");
    SparsityPattern sparsity_pattern_viscosity_matrix;
    sparsity_pattern_viscosity_matrix.init(num_rows);

    scalar *mem_loc;
    matrix3 J = {};

    // Create a temporary lookup for checking if a node is 'pinned' or not.
    // if it is, then only a 1 on the diagonal corresponding to that node should
    // be placed (no off diagonal), effectively taking this node out of the equation
    // and therefore meaning the force on it should always be zero.
    vector<int> is_pinned(node.size(), 0);
    for (int i = 0; i < pinned_nodes_list.size(); ++i) {
        is_pinned[pinned_nodes_list[i]] = 1;
    }
    for(set<int>::iterator it = bsite_pinned_node_list.begin(); it != bsite_pinned_node_list.end(); ++it) {
        is_pinned[*it] = 1;
    }

    for (int n = 0; n < elem.size(); n++) {
        elem[n].calculate_jacobian(J);
        elem[n].calc_shape_function_derivatives_and_volume(J);
        elem[n].create_viscosity_matrix();
        for (int ni = 0; ni < 10; ++ni) {
            for (int nj = 0; nj < 10; ++nj) {
                int ni_index = elem[n].n[ni]->index;
                int nj_index = elem[n].n[nj]->index;
                int ni_row = ni_index * 3;
                int nj_row = nj_index * 3;
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
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

    if (params.calc_stokes == 1) {
        for (int ni = 0; ni < num_nodes; ++ni) {
	    if(is_pinned[ni] == 0) {
                for (int nj = 0; nj < 3; ++nj) {
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
    //throw FFEAException();
    // create a preconditioner for solving in less iterations
    // Create the jacobi preconditioner matrix (diagonal)
    try {
        preconditioner = std::vector<scalar>(num_rows);
    } catch(std::bad_alloc &) {
        throw FFEAException("Failed to allocate 'preconditioner' in NoMassCGSolver\n");
    }

    // create the work vectors necessary for use by the conjugate gradient solver in 'solve'
    try {
        r = std::vector<arr3>(num_nodes);
        p = std::vector<arr3>(num_nodes);
        z = std::vector<arr3>(num_nodes);
        q = std::vector<arr3>(num_nodes);
        f = std::vector<arr3>(num_nodes);
    } catch(std::bad_alloc &) {
        throw FFEAException(" Failed to create the work vectors necessary for NoMassCGSolver\n");
    }
}

/*  */
void NoMassCGSolver::solve(std::vector<arr3> &x) {
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
            return;
        }
	//cout << residual2() << " " << epsilon2 << endl;
        delta_old = delta_new;
        delta_new = parallel_apply_preconditioner();
        vec3_scale_and_add(p, z, (delta_new / delta_old), num_nodes);
    }
    // If desired convergence was not reached in the set number of iterations...
    throw FFEAException("Conjugate gradient solver: Could not converge after %d iterations.\n\tEither epsilon or max_iterations_cg are set too low, or something went wrong with the simulation.\n", i_max);
}

/* */
void NoMassCGSolver::print_matrices(std::vector<arr3> &force) {
    V->print_dense_to_file(force);
}

/* */
scalar NoMassCGSolver::conjugate_gradient_residual_assume_x_zero(std::vector<arr3> &b) {
    scalar delta_new = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) shared(b) reduction(+:delta_new)
#endif
    for (int i = 0; i < num_nodes; i++) {
        r[i][0] = b[i][0];
        r[i][1] = b[i][1];
        r[i][2] = b[i][2];
        f[i][0] = b[i][0];
        f[i][1] = b[i][1];
        f[i][2] = b[i][2];
        b[i][0] = 0;
        b[i][1] = 0;
        b[i][2] = 0;
        z[i][0] = preconditioner[(3 * i)] * r[i][0];
        z[i][1] = preconditioner[(3 * i) + 1] * r[i][1];
        z[i][2] = preconditioner[(3 * i) + 2] * r[i][2];
        p[i][0] = z[i][0];
        p[i][1] = z[i][1];
        p[i][2] = z[i][2];
        delta_new += r[i][0] * z[i][0] + r[i][1] * z[i][1] + r[i][2] * z[i][2];
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
        r2 += r[i][0] * r[i][0] + r[i][1] * r[i][1] + r[i][2] * r[i][2];
        f2 += f[i][0] * f[i][0] + f[i][1] * f[i][1] + f[i][2] * f[i][2];
    }
    if (f2 == 0.0) {
        return 0.0;
    } else {
        return r2 / f2;
    }
}

/* */
scalar NoMassCGSolver::modx(const std::vector<arr3> &x) {
    scalar r2 = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) shared(x) reduction(+:r2)
#endif
    for (int i = 0; i < num_nodes; i++) {
        r2 += x[i][0] * x[i][0] + x[i][1] * x[i][1] + x[i][2] * x[i][2];
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
        pTq += p[i][0] * q[i][0] + p[i][1] * q[i][1] + p[i][2] * q[i][2];
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
        z[i][0] = preconditioner[(3 * i)] * r[i][0];
        z[i][1] = preconditioner[(3 * i) + 1] * r[i][1];
        z[i][2] = preconditioner[(3 * i) + 2] * r[i][2];
        delta_new += r[i][0] * z[i][0] + r[i][1] * z[i][1] + r[i][2] * z[i][2];
    }

    return delta_new;
}

/* */
void NoMassCGSolver::check(const std::vector<arr3> &x) {
    FILE *fout2;
    fout2 = fopen("/localhome/py09bh/output/nomass/cube_viscosity_no_mass.csv", "a");
    int i;
    double temp = 0, temp2 = 0;
    std::vector<arr3> temp_vec = std::vector<arr3>(num_nodes);
    V->apply(x, temp_vec);
    for (i = 0; i < num_nodes; ++i) {
        temp += x[i][0] * temp_vec[i][0] + x[i][1] * temp_vec[i][1] + x[i][2] * temp_vec[i][2];
        temp2 += x[i][0] * f[i][0] + x[i][1] * f[i][1] + x[i][2] * f[i][2];
    }
    fprintf(fout2, "%e,%e\n", temp2, fabs(temp - temp2));
    fclose(fout2);
}
