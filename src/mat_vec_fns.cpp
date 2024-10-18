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

#include "mat_vec_fns.h"

void mat12_apply(matrix12 &A, vector12 &v) {
    int i, j;
    scalar temp_v[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    for (i = 0; i < 12; i++)
        for (j = 0; j < 12; j++)
            temp_v[i] += A[i][j] * v[j];
    for (i = 0; i < 12; i++) v[i] = temp_v[i];
}

void vec3_mat3_mult(arr3 &v, matrix3 &A, arr3 &notv) {
    //int i, j;

    notv[0] = A[0][0]*v[0] + A[1][0]*v[1] + A[2][0]*v[2];
    notv[1] = A[0][1]*v[0] + A[1][1]*v[1] + A[2][1]*v[2];
    notv[2] = A[0][2]*v[0] + A[1][2]*v[1] + A[2][2]*v[2];

}

/* void vec3_vec3_subs(arr3 *u, arr3 *v, arr3 *w) {
    
    w->x = u->x - v->x;
    w->y = u->y - v->y;
    w->z = u->z - v->z;

}

void vec3_vec3_cross(arr3 *u, arr3 *v, arr3 *w) {

	w->x = u->y * v->z - u->z * v->y;
	w->y = u->z * v->x - u->x * v->z;
	w->z = u->x * v->y - u->y * v->x;
}*/ 

/*
 *
 */
void mat3_mult(matrix3 &A, matrix3 &B, matrix3 &result) {
    int i, j, k;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            for (k = 0; k < 3; k++)
                result[i][j] += A[i][k] * B[k][j];
}

/*
 *
 */
void vec12_add(vector12 &A, vector12 &B) {
    for (int i = 0; i < 12; ++i) {
        B[i] = A[i] + B[i];
    }
}

/*
 *
 */
void vec12_scale(vector12 &A, scalar scale) {
    for (int i = 0; i < 12; ++i) {
        A[i] = A[i] * scale;
    }
}

/*
 *
 */
void mat3_mult_transpose(matrix3 &A, matrix3 &B, matrix3 &result) {
    int i, j, k;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            for (k = 0; k < 3; k++)
                result[i][j] += A[i][k] * B[j][k];
}

/*
 *
 */
void mat3_mult_both_transposed(matrix3 &A, matrix3 &B, matrix3 &result) {
    int i, j, k;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            for (k = 0; k < 3; k++)
                result[i][j] += A[k][i] * B[j][k];
}

/*
 *
 */
void mat3_scale(matrix3 &A, scalar s) {
    A[0][0] *= s;
    A[0][1] *= s;
    A[0][2] *= s;
    A[1][0] *= s;
    A[1][1] *= s;
    A[1][2] *= s;
    A[2][0] *= s;
    A[2][1] *= s;
    A[2][2] *= s;
}

/*
 *
 */
scalar mat3_double_contraction_symmetric(matrix3 &A) {
    return A[0][0] * A[0][0] + A[1][1] * A[1][1] + A[2][2] * A[2][2]
            + 2 * (A[0][1] * A[0][1] + A[0][2] * A[0][2] + A[1][2] * A[1][2]);
}

/*
 *
 */
scalar mat3_double_contraction(matrix3 &A) {
    return A[0][0] * A[0][0] + A[1][1] * A[1][1] + A[2][2] * A[2][2]
            + A[0][1] * A[0][1] + A[1][0] * A[1][0] + A[0][2] * A[0][2]
            + A[2][0] * A[2][0] + A[1][2] * A[1][2] + A[2][1] * A[2][1];
}

/*
 *
 */
void mat3_invert(matrix3 &m, matrix3 &m_inv, scalar *det_m) {
    scalar det;

    // Construct the inverse matrix
    m_inv[0][0] = m[2][2] * m[1][1] - m[2][1] * m[1][2];
    m_inv[0][1] = m[2][1] * m[0][2] - m[2][2] * m[0][1];
    m_inv[0][2] = m[1][2] * m[0][1] - m[1][1] * m[0][2];
    m_inv[1][0] = m[2][0] * m[1][2] - m[2][2] * m[1][0];
    m_inv[1][1] = m[2][2] * m[0][0] - m[2][0] * m[0][2];
    m_inv[1][2] = m[1][0] * m[0][2] - m[1][2] * m[0][0];
    m_inv[2][0] = m[2][1] * m[1][0] - m[2][0] * m[1][1];
    m_inv[2][1] = m[2][0] * m[0][1] - m[2][1] * m[0][0];
    m_inv[2][2] = m[1][1] * m[0][0] - m[1][0] * m[0][1];

    // calc determinant
    det = m[0][0] * m_inv[0][0] + m[1][0] * m_inv[0][1] + m[2][0] * m_inv[0][2];
    *det_m = det;

    // divide by determinant
    det = 1.0 / det;
    m_inv[0][0] *= det;
    m_inv[0][1] *= det;
    m_inv[0][2] *= det;
    m_inv[1][0] *= det;
    m_inv[1][1] *= det;
    m_inv[1][2] *= det;
    m_inv[2][0] *= det;
    m_inv[2][1] *= det;
    m_inv[2][2] *= det;
}

/*
 *
 */
void mat12_set_zero(matrix12 &A) {
    int i, j;
    for (i = 0; i < 12; i++)
        for (j = 0; j < 12; j++)
            A[i][j] = 0;
}

/*
 *
 */
void mat3_set_zero(matrix3 &A) {
    int i, j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            A[i][j] = 0;
}

/*
 *
 */
void mat3_set_identity(matrix3 &A) {

    mat3_set_zero(A);
    int i;
    for (i = 0; i < 3; i++)
        A[i][i] = 1;
}

/*
 *
 */
void mat4_set_zero(matrix4 &A) {
    int i, j;
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            A[i][j] = 0;
}

/*
 *
 */
void arr3_set_zero(arr3 &v) {
   
    v[0] = 0;
    v[1] = 0;
    v[2] = 0;

}

/*
 *
void vec3_scale(std::vector<arr3> &v, scalar scale) {
	v->x *= scale;
	v->y *= scale;
	v->z *= scale;
}

void vec3_scale2(std::vector<arr3> &v1, std::vector<arr3> &v2, scalar scale) {
	v2->x = scale*v1->x;
	v2->y = scale*v1->y;
	v2->z = scale*v1->z;
}
 */


/*
 *
 */
void vec3_add_to_scaled(std::vector<arr3> &v1, std::vector<arr3> &v2, scalar a, int vec_size) {
    int i;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i) shared(v1, v2, a, vec_size)
#endif
    for (i = 0; i < vec_size; i++) {
        v1[i][0] += a * v2[i][0];
        v1[i][1] += a * v2[i][1];
        v1[i][2] += a * v2[i][2];
    }
}

/*
 *
 */
void vec3_scale_and_add(std::vector<arr3> &v1, std::vector<arr3> &v2, scalar a, int vec_size) {
    int i;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i) shared(v1, v2, a, vec_size)
#endif
    for (i = 0; i < vec_size; i++) {
        v1[i][0] = a * v1[i][0] + v2[i][0];
        v1[i][1] = a * v1[i][1] + v2[i][1];
        v1[i][2] = a * v1[i][2] + v2[i][2];
    }
}

/*
 *
 */
void vec12_set_zero(vector12 &v) {
    int i;
    for (i = 0; i < 12; i++)
        v[i] = 0;
}

/*
 *
 */
void print_matrix3(matrix3 &m) {
    int i;
    for (i = 0; i < 3; i++)
        printf("%e %e %e\n", m[i][0], m[i][1], m[i][2]);
}

/*
 *
 */
void print_matrix4(matrix4 &m) {
    int i;
    for (i = 0; i < 4; i++)
        printf("%e %e %e %e\n", m[i][0], m[i][1], m[i][2], m[i][3]);
}

/*
 *
 */
void print_matrix12(matrix12 &m) {
    int i, j;
    for (i = 0; i < 12; i++) {
        for (j = 0; j < 12; j++)
            printf("%e ", m[i][j]);
        printf("\n");
    }
}

/*
 * 
 */
void print_vector12(vector12 &v) {
    int i;
    for (i = 0; i < 12; i++)
        printf("%e\n", v[i]);
}

void print_arr3(arr3 &v) {
    printf("%e %e %e\n", v[0], v[1], v[2]);
}

/*scalar mag(arr3 &v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

std::array<scalar,3> normalise(arr3 &v) {
    scalar magnitude;
    arr3 norm;
    magnitude = mag(v);
    if(magnitude == 0.0) {
	throw -1;
    }
    norm.x = v.x / magnitude;
    norm.y = v.y / magnitude;
    norm.z = v.z / magnitude;
    return norm.data;
} */ 
