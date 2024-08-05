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
 *      BEM_Poisson_Boltzmann.h
 *	Author: Robin Richardson, University of Leeds
 *	Email: pyrar@leeds.ac.uk
 */

#define _USE_MATH_DEFINES
#include <cmath>

#include "FFEA_return_codes.h"
#include "BEM_Poisson_Boltzmann.h"
#include "NearestNeighbourLinkedListCube.h"
#include "SparseMatrixUnknownPattern.h"

#include "GaussianQuadrature_tri.h"
#include "GaussianQuadrature_1d.h"

BEM_Poisson_Boltzmann::BEM_Poisson_Boltzmann() {
    num_faces = 0;
    lookup = NULL;
    kappa = 0;
    mat_C = NULL;
    mat_D = NULL;
}

BEM_Poisson_Boltzmann::~BEM_Poisson_Boltzmann() {
    num_faces = 0;
    lookup = NULL;
    kappa = 0;

    delete mat_C;
    delete mat_D;
    mat_C = NULL;
    mat_D = NULL;
}

int BEM_Poisson_Boltzmann::init(NearestNeighbourLinkedListCube *lookup) {
    this->lookup = lookup;
    this->num_faces = lookup->get_pool_size();

    //			int n2 = num_faces * num_faces;
    //			mat_C = new scalar[n2];
    //			mat_D = new scalar[n2];

    /* Create and initialise our sparse matrices */
    mat_C = new(std::nothrow) SparseMatrixUnknownPattern();
    if (mat_C == NULL) FFEA_ERROR_MESSG("Could not allocate C matrix\n");
    if (mat_C->init(num_faces, 100) == FFEA_ERROR) {
        FFEA_ERROR_MESSG("Could not allocate memory for C matrix\n")
    }

    mat_D = new(std::nothrow) SparseMatrixUnknownPattern();
    if (mat_D == NULL) FFEA_ERROR_MESSG("Could not allocate D matrix\n");
    if (mat_D->init(num_faces, 100) == FFEA_ERROR) {
        FFEA_ERROR_MESSG("Could not allocate memory for D matrix\n")
    }

    return FFEA_OK;
}

/*
 * Sets the inverse debye screening length for the system
 */
void BEM_Poisson_Boltzmann::set_kappa(scalar kappa) {
    this->kappa = kappa;
}

/*
 *
 */
void BEM_Poisson_Boltzmann::build_BEM_matrices() {
    int i;
    LinkedListNode<Face> *l_i;
    Face *f;
    vector3 gqp[4];

    /* Clear the C and D matrices */
    mat_C->zero();
    mat_D->zero();

    /* For each face, calculate the interaction with all other relevant faces and add the contribution to mat_C and mat_D */
    for (i = 0; i < num_faces; i++) {

        // get the ith face
        l_i = lookup->get_from_pool(i);
        f = l_i->obj;

        // Create matrix C diagonal (self term) for constant element case
        mat_C->set_diagonal_element(i, -.5 * 4.0 * M_PI);

        // Create matrix D diagonal (self term) for constant element case
        mat_D->set_diagonal_element(i, -(self_term(&f->centroid, &f->n[1]->pos, &f->n[2]->pos, 6) +
                self_term(&f->centroid, &f->n[2]->pos, &f->n[0]->pos, 6) +
                self_term(&f->centroid, &f->n[0]->pos, &f->n[1]->pos, 6)));


        // calculate the 4 gauss points for this face
        f->barycentric_calc_point(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, &gqp[0]);
        f->barycentric_calc_point(.6, .2, .2, &gqp[1]);
        f->barycentric_calc_point(.2, .6, .2, &gqp[2]);
        f->barycentric_calc_point(.2, .2, .6, &gqp[3]);

        // Perform the necessary integrals for the face f with all faces in its own cell (in the lookup grid)
        perform_integrals_for_lookup_cell_self(l_i, gqp);

        // Perform integrals for the 9 adjacent cells above and below the plane containg f's cell
        for (int x = -1; x < 2; x++)
            for (int z = -1; z < 2; z++) {
                perform_integrals_for_lookup_cell_relative(l_i, gqp, x, -1, z);
                perform_integrals_for_lookup_cell_relative(l_i, gqp, x, +1, z);
            }

        // Perform integrals for the remaining 8 adjacent cells, lying in the y-plane of f's cell
        for (int x = -1; x < 2; x++) {
            perform_integrals_for_lookup_cell_relative(l_i, gqp, x, 0, -1);
            perform_integrals_for_lookup_cell_relative(l_i, gqp, x, 0, +1);
        }
        perform_integrals_for_lookup_cell_relative(l_i, gqp, -1, 0, 0);
        perform_integrals_for_lookup_cell_relative(l_i, gqp, +1, 0, 0);
    }
}

void BEM_Poisson_Boltzmann::perform_integrals_for_lookup_cell_self(LinkedListNode<Face> *l_i, vector3 gqp[4]) {
    scalar mat_C_contribution, mat_D_contribution;

    // Get the top of the stack in the cell in which face f lies
    LinkedListNode<Face> *l_j = lookup->get_top_of_stack(l_i->x, l_i->y, l_i->z);
    while (l_j != NULL) {
        if (l_j->index != l_i->index) {
            gauss_quadrature_4_point(gqp,
                    &(l_j->obj->centroid),
                    &mat_D_contribution,
                    &mat_C_contribution,
                    l_i->obj);

            mat_C->add_off_diagonal_element(l_i->index, l_j->index, mat_C_contribution);
            mat_D->add_off_diagonal_element(l_i->index, l_j->index, mat_D_contribution);
        }
        l_j = l_j->next;
    }
}

void BEM_Poisson_Boltzmann::perform_integrals_for_lookup_cell_relative(LinkedListNode<Face> *l_i, vector3 gqp[4], int dx, int dy, int dz) {
    scalar mat_C_contribution, mat_D_contribution;

    LinkedListNode<Face> *l_j = lookup->get_top_of_stack(l_i->x + dx, l_i->y + dy, l_i->z + dz);
    while (l_j != NULL) {
        gauss_quadrature_4_point(gqp,
                &(l_j->obj->centroid),
                &mat_D_contribution,
                &mat_C_contribution,
                l_i->obj);

        mat_C->add_off_diagonal_element(l_i->index, l_j->index, mat_C_contribution);
        mat_D->add_off_diagonal_element(l_i->index, l_j->index, mat_D_contribution);

        l_j = l_j->next;
    }
}

void BEM_Poisson_Boltzmann::print_matrices() {
    printf("\n\nmat_C:\n\n");
    mat_C->print();

    printf("\n\nmat_D:\n\n");
    mat_D->print();
}

SparseMatrixUnknownPattern * BEM_Poisson_Boltzmann::get_C() {
    return mat_C;
}

SparseMatrixUnknownPattern * BEM_Poisson_Boltzmann::get_D() {
    return mat_D;
}

/* Returns the value of the fundamental solution u multiplied by 4*pi */
scalar BEM_Poisson_Boltzmann::u_4pi(scalar r) {
    return exp(-kappa * r) / r;
}

/* Returns the radial component of grad of u multiplied by 4*pi (all other components are zero) */
scalar BEM_Poisson_Boltzmann::grad_u_4pi(scalar r, scalar r2) {
    return -exp(-kappa * r) * (1.0 / r2 + kappa / r);
}

/*
                scalar BEM_Poisson_Boltzmann::screened_R_theta(scalar r_perp_mag, scalar half_theta_max, scalar theta_bar, scalar xi)
                {
                        return exp( (-kappa * r_perp_mag) / cos(half_theta_max * (xi + 1) - theta_bar) );
                }
 */

void BEM_Poisson_Boltzmann::gauss_quadrature_4_point(vector3 gqp[4], vector3 *p, scalar *int_u, scalar *int_du, Face *f) {
    scalar r_hat_dot_n[4], r2[4], r_mag[4];
    vector3 r;

    int i;
    for (i = 0; i < 4; i++) {
        r.x = gqp[i].x - p->x;
        r.y = gqp[i].y - p->y;
        r.z = gqp[i].z - p->z;
        r2[i] = r.x * r.x + r.y * r.y + r.z * r.z;
        r_mag[i] = sqrt(r2[i]);
        r_hat_dot_n[i] = (r.x * f->normal.x + r.y * f->normal.y + r.z * f->normal.z) / r_mag[i];
    }

    *int_u = -f->area * (-.5625 * u_4pi(r_mag[0]) + 0.5208333333 * (u_4pi(r_mag[1]) + u_4pi(r_mag[2]) + u_4pi(r_mag[3])));
    *int_du = f->area * (-.5625 * grad_u_4pi(r_mag[0], r2[0]) * r_hat_dot_n[0] + 0.5208333333 * (
            grad_u_4pi(r_mag[1], r2[1]) * r_hat_dot_n[1] +
            grad_u_4pi(r_mag[2], r2[2]) * r_hat_dot_n[2] +
            grad_u_4pi(r_mag[3], r2[3]) * r_hat_dot_n[3]));
}

scalar BEM_Poisson_Boltzmann::self_term(vector3 *n0, vector3 *n1, vector3 *n2, int precision) {
    // Calculate the triangle side vectors, r_01, r_02 and r_12
    vector3 r_01, r_02, r_12; // r_12;
    r_01.assign ( n1->x - n0->x, n1->y - n0->y, n1->z - n0->z );
    r_02.assign ( n2->x - n0->x, n2->y - n0->y, n2->z - n0->z );
    r_12.assign ( n2->x - n1->x, n2->y - n1->y, n2->z - n1->z );

    // Get the lengths of these vectors
    scalar r_01_sq = r_01.x * r_01.x + r_01.y * r_01.y + r_01.z * r_01.z,
            r_02_sq = r_02.x * r_02.x + r_02.y * r_02.y + r_02.z * r_02.z,
            r_12_sq = r_12.x * r_12.x + r_12.y * r_12.y + r_12.z * r_12.z;

    // Calculate some dot products
    scalar r_01_dot_r_02 = r_01.x * r_02.x + r_01.y * r_02.y + r_01.z * r_02.z,
            r_01_dot_r_12 = r_01.x * r_12.x + r_01.y * r_12.y + r_01.z * r_12.z;

    // Get the size of the angle through which we need to integrate
    scalar theta_max = acos(r_01_dot_r_02 / sqrt(r_01_sq * r_02_sq));

    // Get the vector between node 0 and the point perpendicularly opposite (on line node 1 to 2)
    vector3 r_0_perp;
    r_0_perp.assign ( r_01.x - r_01_dot_r_12 * r_12.x / r_12_sq,
        r_01.y - r_01_dot_r_12 * r_12.y / r_12_sq,
        r_01.z - r_01_dot_r_12 * r_12.z / r_12_sq );

    // Get the perpendicular distance between node 0 and line 1 to 2
    scalar L_perp = sqrt(r_0_perp.x * r_0_perp.x + r_0_perp.y * r_0_perp.y + r_0_perp.z * r_0_perp.z);

    // Get the angle between side 01 and the perpedicular line between 0 and side 12.
    scalar theta_star = acos(L_perp / sqrt(r_01_sq));

    return integrate_function_1d_tri(theta_max, L_perp, theta_star, precision);
}

scalar BEM_Poisson_Boltzmann::f_1d(scalar r) {
    return exp(-kappa * r);
}

scalar BEM_Poisson_Boltzmann::f_3d(vector3 *p, vector3 *q) {
    scalar dx = p->x - q->x, dy = p->y - q->y, dz = p->z - q->z;
    scalar r = sqrt(dx * dx + dy * dy + dz * dz);
    return exp(-kappa * r) / r;
}

