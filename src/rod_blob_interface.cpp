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
 */
 
// Note: the code herein will often refer to the attachment element as
// the 'attachment_node' and the attachmen node as the
// 'attachment_node_pos'. This is a hangover from an earlier stage of
// development, and will eventually be fixed, but please keep it in mind.
// Comments will clarify where appropriate.

#include "rod_blob_interface.h" 

namespace rod{

// interface math

/**
 Given the three nodes in a triangle (or the face of a tetrahedron, if you prefer), get the direction vector of the normal to to that face. 
*/ 
void get_tri_norm(float node0[3], float node1[3], float node2[3], OUT float tri_norm[3]){
    float ax = node1[0] - node0[0];
	float ay = node1[1] - node0[1];
	float az = node1[2] - node0[2];
	float bx = node2[0] - node1[0];
	float by = node2[1] - node1[1];
	float bz = node2[2] - node1[2];
    tri_norm[0] = az*by - ay*bz;
    tri_norm[1] = ax*bz - az*bx;
    tri_norm[2] = ay*bx - ax*by;
}

/**
 For an array of 9 tet_nodes (defined in  object defined in mesh_node.cpp/mesh_node.h) populate a 1-D array with the 3x3 Jacobian of the shape functions. 
*/
void get_jacobian(mesh_node **tet_nodes, OUT float J[9]){
    J[0] = tet_nodes[1]->pos[0] - tet_nodes[0]->pos[0];
    J[1] = tet_nodes[1]->pos[1] - tet_nodes[0]->pos[1];
    J[2] = tet_nodes[1]->pos[2] - tet_nodes[0]->pos[2];
    
    J[3] = tet_nodes[2]->pos[0] - tet_nodes[0]->pos[0];
    J[4] = tet_nodes[2]->pos[1] - tet_nodes[0]->pos[1];
    J[5] = tet_nodes[2]->pos[2] - tet_nodes[0]->pos[2];
    
    J[6] = tet_nodes[3]->pos[0] - tet_nodes[0]->pos[0];
    J[7] = tet_nodes[3]->pos[1] - tet_nodes[0]->pos[1];
    J[8] = tet_nodes[3]->pos[2] - tet_nodes[0]->pos[2];
    
    
    //J[0] = tet_nodes[1]->pos[0] - tet_nodes[0]->pos[0];
    //J[3] = tet_nodes[1]->pos[1] - tet_nodes[0]->pos[1];
    //J[6] = tet_nodes[1]->pos[2] - tet_nodes[0]->pos[2];
    
    //J[1] = tet_nodes[2]->pos[0] - tet_nodes[0]->pos[0];
    //J[4] = tet_nodes[2]->pos[1] - tet_nodes[0]->pos[1];
    //J[7] = tet_nodes[2]->pos[2] - tet_nodes[0]->pos[2];
    
    //J[2] = tet_nodes[3]->pos[0] - tet_nodes[0]->pos[0];
    //J[5] = tet_nodes[3]->pos[1] - tet_nodes[0]->pos[1];
    //J[8] = tet_nodes[3]->pos[2] - tet_nodes[0]->pos[2];
    
}

/**
 Find the inverse of a 3x3 matrix (m), and populate the 1-D array m_inv. This is identical to the matrix inversion in mat_vec_fns, only it operates on 1-D arrays of floats instead of 2-D ones of scalars.
*/ 
void float_3x3_invert(float m[9], OUT float m_inv[9]) {
    scalar det;

    // Construct the inverse matrix
    m_inv[0] = m[8] * m[4] - m[7] * m[5];
    m_inv[1] = m[7] * m[2] - m[8] * m[1];
    m_inv[2] = m[5] * m[1] - m[4] * m[2];
    m_inv[3] = m[6] * m[5] - m[8] * m[3];
    m_inv[4] = m[8] * m[0] - m[6] * m[2];
    m_inv[5] = m[3] * m[2] - m[5] * m[0];
    m_inv[6] = m[7] * m[3] - m[6] * m[4];
    m_inv[7] = m[6] * m[1] - m[7] * m[0];
    m_inv[8] = m[4] * m[0] - m[3] * m[1];

    // calc determinant
    det = m[0] * m_inv[0] + m[3] * m_inv[1] + m[6] * m_inv[2];

    // divide by determinant
    det = 1.0 / det;
    m_inv[0] *= det;
    m_inv[1] *= det;
    m_inv[2] *= det;
    m_inv[3] *= det;
    m_inv[4] *= det;
    m_inv[5] *= det;
    m_inv[6] *= det;
    m_inv[7] *= det;
    m_inv[8] *= det;
}

/**
 Multiply the 3x3 matrices A and B. This one is unrolled. Store the result in result. A, B and result are 1-D arrays of length 9 representing 3x3 matrices. Again, this does the same thing as the mat_vec_fn, but with 1-D arrays of floats.
*/ 
void float_3x3_mult_unrolled(float A[9], float B[9], OUT float result[9]){
    result[0] += A[0] * B[0];
    result[0] += A[3] * B[1];
    result[0] += A[6] * B[2];
    result[1] += A[0] * B[3];
    result[1] += A[3] * B[4];
    result[1] += A[6] * B[5];
    result[2] += A[0] * B[6];
    result[2] += A[3] * B[7];
    result[2] += A[6] * B[8];
    result[3] += A[1] * B[0];
    result[3] += A[4] * B[1];
    result[3] += A[7] * B[2];
    result[4] += A[1] * B[3];
    result[4] += A[4] * B[4];
    result[4] += A[7] * B[5];
    result[5] += A[1] * B[6];
    result[5] += A[4] * B[7];
    result[5] += A[7] * B[8];
    result[6] += A[2] * B[0];
    result[6] += A[5] * B[1];
    result[6] += A[8] * B[2];
    result[7] += A[2] * B[3];
    result[7] += A[5] * B[4];
    result[7] += A[8] * B[5];
    result[8] += A[2] * B[6];
    result[8] += A[5] * B[7];
    result[8] += A[8] * B[8];
}

/**
 Transpose a 3x3 matrix, represented as a 1-D array of floats. Populate the array 'transposed'.
*/
void transpose_3x3(float in[9], OUT float transposed[9]){
    transposed[0] = in[0];
    transposed[1] = in[3];
    transposed[2] = in[6];
    transposed[3] = in[1];
    transposed[4] = in[4];
    transposed[5] = in[7];
    transposed[6] = in[2];
    transposed[7] = in[5];
    transposed[8] = in[8];
}

// note: this is mostly the same as the one in tetra_element_linear, but that one is a method that changes
// the member variables of an object. i can't use that one here because i specifically don't WANT to change
// those variables. so thanks bjarne
/**
 Get the 3x3 gradient deformation matrix of a tetrahedron. The tetrahedron is given as an array of the mesh_node object defined in mesh_node.cpp. Two arrays are needed, one representing the tetrahedron before it was deformed, the other one representing the tetrahedron after deformation.
 \f[ \mathbf{F_e} = (\mathbf{J'} \mathbf{J}^{-1})^T \f]
 Where \f$ \mathbf{F} \f$ is the gradient deformation matrix, \f$ \mathbf{J'} \f$ is the Jacobian of the new deformed tetrahedron, and \f$ \mathbf{J} \f$ is the Jacobian of the undeformed one.
 Note: this is mostly the same as the one in tetra_element_linear.cpp, but that one is a member function that operates on member variables of tetra_element_linear, whereas this one is a pure function (and uses the same floats and 1-d arrays as the rest of the rod stuff).
*/
void get_gradient_deformation(float J_inv_0[9], mesh_node**nodes_curr, OUT float transposed_gradient_deformation_3x3[9]){
    if(dbg_print){std::cout << " Gradient deformation is occuring.\n";}
    if(dbg_print){std::cout << "  Getting jacobian...\n";}
    //float J_before[9];
    float J_after[9];
    //get_jacobian(nodes_before, J_before);
    get_jacobian(nodes_curr, J_after);
    //float before_inverse[9];
    //float_3x3_invert(J_before, before_inverse);
    float gradient_deformation_3x3[9] = {0,0,0,0,0,0,0,0,0};
    float_3x3_mult_unrolled(J_after, J_inv_0, gradient_deformation_3x3);
    transpose_3x3(gradient_deformation_3x3, transposed_gradient_deformation_3x3);
    
//    print_array("  Tet before node 0", nodes_before[0]->pos.data, 3);
//    print_array("  Tet before node 1", nodes_before[1]->pos.data, 3);
//    print_array("  Tet before node 2", nodes_before[2]->pos.data, 3);
//    print_array("  Tet before node 3", nodes_before[3]->pos.data, 3);
    print_array("  Tet after node 0", nodes_curr[0]->pos.data, 3);
    print_array("  Tet after node 1", nodes_curr[1]->pos.data, 3);
    print_array("  Tet after node 2", nodes_curr[2]->pos.data, 3);
    print_array("  Tet after node 3", nodes_curr[3]->pos.data, 3);
    
    print_array("  J after", J_after, 9);
    print_array("  before inverse", J_inv_0, 9);
    print_array("  gradient_deformation_3x3", gradient_deformation_3x3, 9);
    print_array("  transposed_gradient_deformation_3x3", transposed_gradient_deformation_3x3, 9);
}

/**
 QR decompose\factorise a 3x3 matrix (given here as a 1-d array of floats) using the Gram-Schmidt process. Like the matrix multiplication above, this is again unrolled. This populates arrays for both Q and R, although only R is used.
 We wish to decompose the 3x3 matrix A ([a_1 ... a_n]) into two 3x3 matrices, Q and R.
 The elements e_n of our Q matrix are given by the following recursive definition:
 
 \f[ u_1 = a_1, ~ e_1 = \frac{u1}{||u_1||} \f]
 \f[ u_{k+1} = a_{k+1} - (a_{k+1} \cdot e_1)e1 - ... - (a_{k+1} \cdot e_k)e_k, ~ e_1 = \frac{u1}{||u_1||}, e_{k+1} = \frac{u_{k+1}}{||u_{k+1}||} \f]
 
 \f[A = [ a_1 | a_2 | \dots | a_n] = [ e_1 | e_2 | \dots | e_n]
    \begin{bmatrix} 
    a_{1}\cdot e_1 & a_2\cdot e_1 & \dots & a_n \cdot e_1 \\
    0 & a_2 \cdot e_2 & \dots & a_n \cdot e_2 \\ 
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \dots & a_n \cdot e_n 
    \end{bmatrix}

  = QR \f]
*/ 
void QR_decompose_gram_schmidt(float matrix_3x3[9], OUT float Q[9], float R[9]){
    float a1[3] = {matrix_3x3[0], matrix_3x3[3], matrix_3x3[6]};
    float a2[3] = {matrix_3x3[1], matrix_3x3[4], matrix_3x3[7]};
    float a3[3] = {matrix_3x3[2], matrix_3x3[5], matrix_3x3[8]};
    float u1[3] = {a1[0], a1[1], a1[2]};
    float e1[3];
    normalize(u1, e1);
    float u2[3];
    float a2_dot_e1 = a2[0]*e1[0] + a2[1]*e1[1] + a2[2]*e1[2];
    vec3d(n){ u2[n] = a2[n] - a2_dot_e1*e1[n];}
    float e2[3];
    normalize(u2, e2);
    float u3[3];
    float a3_dot_e1 = a3[0]*e1[0] + a3[1]*e1[1] + a3[2]*e1[2];
    float a3_dot_e2 = a3[0]*e2[0] + a3[1]*e2[1] + a3[2]*e2[2];
    vec3d(n) u3[n] = a3[n] - a3_dot_e1*e1[n] - a3_dot_e2*e2[n];
    float e3[3];
    normalize(u3, e3);
    
    Q[0] = e1[0]; Q[1] = e1[1]; Q[2] = e1[2]; Q[3] = e2[0]; Q[4] = e2[1]; Q[5] = e2[2]; Q[6] = e3[0]; Q[7] = e3[1]; Q[8] = e3[2];
    R[0] = a1[0]*e1[0] + a1[1]*e1[1] + a1[2]*e1[2];
    R[1] = a2[0]*e1[0] + a2[1]*e1[1] + a2[2]*e1[2];
    R[2] = a3[0]*e1[0] + a3[1]*e1[1] + a3[2]*e1[2];
    R[3] = 0;
    R[4] = a2[0]*e2[0] + a2[1]*e2[1] + a2[2]*e2[2];
    R[5] = a3[0]*e2[0] + a3[1]*e2[1] + a3[2]*e2[2];
    R[6] = 0;
    R[7] = 0;
    R[8] = a3[0]*e3[0] + a3[1]*e3[1] + a3[2]*e3[2];
    
    //float Q[9] = { e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], e3[0], e3[1], e3[2] }
    //float R[9] = { 
    //    a1[0](e1[0] + a1[1](e1[1] + a1[2](e1[2], a2[0](e1[0] + a2[1](e1[1] + a2[2](e1[2], a3[0](e1[0] + a3[1](e1[1] + a3[2](e1[2],
    //    0,  a2[0](e2[0] + a2[1](e2[1] + a2[2](e2[2], a3[0](e2[0] + a3[1](e2[1] + a3[2](e2[2]
    //    0, 0, a3[0](e3[0] + a3[1](e3[1] + a3[2](e3[2]
}

// X1Y2Z3 :)
/**
 Construct an euler rotation matrix of the X1Y2Z3 type. Parameters, a, b and c, the angles of rotation. Populates the array 'rotmat'.
*/
void construct_euler_rotation_matrix(float a, float b, float g, float rotmat[9]){
    float s1 = std::sin(a);
    float s2 = std::sin(b);
    float s3 = std::sin(g);
    float c1 = std::cos(a);
    float c2 = std::cos(b);
    float c3 = std::cos(g);
    rotmat[0] = c2*c3;
    rotmat[1] = -c2*s3;
    rotmat[2] = s2;
    rotmat[3] = c1*s3 + c3*s1*s2;
    rotmat[4] = c1*c3 - s1*s2*s3;
    rotmat[5] = -c2*s1;
    rotmat[6] = s1*s3 - c1*c3*s2;
    rotmat[7] = c3*s1 + c1*s2*s3;
    rotmat[8] = c1*c2;
}

/**
 Rotate a tetrahedron (given by an array of mesh_node objects) by a 3x3 rotation matrix (given as a 1-d array of floats) about its centroid. Populates a new mesh_node object. Note that this function is not used in the main simulation loop, it is only used in the rod_blob_interface unit test.
*/
void rotate_tet(float rotmat[9], mesh_node **nodes, OUT mesh_node **rotated_nodes){
    
    //print_array("connected tetrahedron node 0", nodes[0]->pos.data, 3);
    //print_array("connected tetrahedron node 1", nodes[1]->pos.data, 3);
    //print_array("connected tetrahedron node 2", nodes[2]->pos.data, 3);
    //print_array("connected tetrahedron node 3", nodes[3]->pos.data, 3);
    
    // get centroid
    float tet_centroid[3] = {0,0,0};
    for (int i=0; i<4; i++){
        vec3d(n) {tet_centroid[n] += nodes[i]->pos.data[n];}
        
    }
    
    vec3d(n){ tet_centroid[n] /= 4.0; }
    //print_array("tet_centroid", tet_centroid, 3);
    
    // move to centre
    for (int i=0; i<4; i++){
        vec3d(n){ nodes[i]->pos[n] -= tet_centroid[n]; }
    }
    
    float posdata_tofloat[3];
    float posdata_rotated[3];

    // rotate
    for (int i=0; i<4; i++){
    //    print_array("posdata_toflnfr", nodes[i]->pos.data, 3);
        vec3d(n){posdata_tofloat[n] = (float)nodes[i]->pos[n];}
        apply_rotation_matrix(posdata_tofloat, rotmat, posdata_rotated);
        vec3d(n){rotated_nodes[i]->pos.data[n] = posdata_rotated[n];}
    //    print_array("rotated_nodes[i]->pos.data[n]", rotated_nodes[i]->pos.data, 3);
    }
    
    // move back
    for (int i=0; i<4; i++){
        vec3d(n){ nodes[i]->pos[n] += tet_centroid[n]; }
        vec3d(n){ rotated_nodes[i]->pos[n] += tet_centroid[n]; }
    }

}
/**
 This is a candidate for removal. It will convert a rotation matrix to
 a set of euler angles. I thought I needed it due to the way the blob
 interface is written, turns out I don't. 
*/
void get_euler_angles(float rm[9], OUT float euler_angles[3]){ // euler angles assumed to be [x, y, z]
    // credit: greg slabaugh: http://www.gregslabaugh.net/publications/euler.pdf
    float theta_1;
    float psi_1;
    float phi_1;
    if ((rm[6] != 1) || (rm[6] != -1)){
        theta_1 = -asin(rm[6]);
        //float theta_2 = M_PI - theta_1;
        psi_1 = atan2(rm[7]/cos(theta_1), rm[8]/cos(theta_1) );
        //float psi_2 = atan2(rm[7]/cos(theta_2), rm[8]/cos(theta_2) );
        phi_1 = atan2(rm[3]/cos(theta_1), rm[0]/cos(theta_1) );
        //float phi_2 = atan2(rm[3]/cos(theta_2), rm[0]/cos(theta_2) );
    }
    else{
        std::cout << "Gimbal lock! Nice!\n";
        float phi_1 = 0;
        if (rm[6] == -1){
            theta_1 = M_PI/2;
            psi_1 = phi_1 + atan2(rm[1], rm[2]);
        }
        else{
            theta_1 = -M_PI/2;
            psi_1 = -phi_1 + atan2(-rm[1], -rm[2]);
        }
    }
    
    euler_angles[x] = phi_1;
    euler_angles[y] = theta_1;
    euler_angles[z] = psi_1;
}

/**
 This one should be pretty self-explanatory, it does the inverse of the
 previous function. Given a set of euler angles, it produces a rotation
 matrix.
*/
void get_rotation_matrix_from_euler(float euler_angles[3], OUT float rm[9]){
    rm[0] = cos(euler_angles[1]) * cos(euler_angles[2]);
    rm[1] = sin(euler_angles[0]) * sin(euler_angles[1]) * cos(euler_angles[2]) - cos(euler_angles[0]) * sin(euler_angles[2]);
    rm[2] = cos(euler_angles[0]) * sin(euler_angles[1]) * cos(euler_angles[2]) + sin(euler_angles[0]) * sin(euler_angles[2]);
    rm[3] = cos(euler_angles[1]) * sin(euler_angles[2]);
    rm[4] = sin(euler_angles[0]) * sin(euler_angles[1]) * sin(euler_angles[2]) + cos(euler_angles[0]) * cos(euler_angles[2]);
    rm[5] = cos(euler_angles[0]) * sin(euler_angles[1]) * sin(euler_angles[2]) - sin(euler_angles[0]) * cos(euler_angles[2]);
    rm[6] = -1 * sin(euler_angles[1]);
    rm[7] = sin(euler_angles[0]) * cos(euler_angles[1]);
    rm[8] = cos(euler_angles[0]) * cos(euler_angles[1]);

}

/**
 Check that all of the elements in two arrays are equal.
 Generally it's not good practice to do comparisons like this based on
 floating point numbers, but this is the only way I have to check for
 equality between face nodes and element nodes (and hence get the indices
 of the face nodes) because whoever programmed that part of the original
 FFEA data structure apparently didn't think that information was
 important.
*/
bool array_equal(float arr1[3], float arr2[3]){
    bool arr_equal = true;
    for (int i=0; i<3; i++){
        if(arr1[i] != arr2[i]){ arr_equal = false; }
    }
    return arr_equal;
}

/**
 Utility function to check if a large [4][3] array contains (non-sequential)
 a smaller [3][3] array. The reason for this code is checking if
 face nodes are contained within the nodes for an element.
 ParamsL:
 * large_arr - the positional data for the nodes of an element.
 * small_arr - the positional data for the nodes of a face.
*/
bool array_contains(float large_arr[4][3], float small_arr[3][3]){
    bool arr_equal = false;
    int num_equal=0;
    for (int i=0; i<4; i++){
        for (int j=0; j<3; j++){
            if (array_equal(large_arr[i], small_arr[j])){
                num_equal+=1;
            }
        }
    }
    if (num_equal == 3){
        arr_equal = true;
    }
    return arr_equal;
}

/**
 Utility function to check if a mesh node is null. 
*/
void mesh_node_null_check(mesh_node* node, std::string location){
    if (!node) { // NULL check
        printf("null node found\n");
        std::cout << "Null node found at " << location << "\n";
    }
}

/**
 This will scale the attachment element such that, if the end element gets
 longer, it gets shorter and vice-versa. The short reason for this is
 so that the attachment element respects the node weighting feature of
 the mutual material axis bend energy calculation, bearing in mind that
 shorter elements have the mutual elements weighted more heavily toward
 them.
 Params:
 * attachment_node, the attachment element
 * end_node, the end element
 * equilibrium versions of those two
 Returns:
 * scaled_attachment_node - the scaled version of the attachment element
 * scaled_attachment_node_equil - the same, at equilibrium
*/
void rescale_attachment_node(float attachment_node[3], float end_node[3], float attachment_node_equil[3], float end_node_equil[3], OUT float scaled_attachment_node[3], float scaled_attachment_node_equil[3]){
    if(dbg_print){std::cout << "      Rescaling attachment node...\n";}
    print_array("      attachment node      ", attachment_node, 3);
    print_array("      attachment node equil", attachment_node_equil, 3);    
    print_array("      end_node             ", end_node, 3);
    print_array("      end_node_equil       ", end_node_equil, 3);
    
    float desired_length = absolute(end_node_equil)*2;
    float attachment_node_scale = desired_length - absolute(end_node);
    float attachment_node_equil_scale = absolute(end_node_equil);
    vec3d(n){scaled_attachment_node[n] = attachment_node[n]*attachment_node_scale;}
    vec3d(n){scaled_attachment_node_equil[n] = attachment_node_equil[n]*attachment_node_equil_scale;}
    if(dbg_print){std::cout << "      desired length       " << desired_length << "\n";}
    if(dbg_print){std::cout << "      attachment_node_scale " << attachment_node_scale << "\n";}
    if(dbg_print){std::cout << "      attachment_node_equil_scale " << attachment_node_equil_scale << "\n";}
    print_array("      scaled_attachment_node", scaled_attachment_node, 3);
    print_array("      scaled_attachment_node_equil", scaled_attachment_node_equil, 3);
}

/**
 Inverts a 3x3 matrix. Yeah, nothing clever going on here.
 Params: m - a matrix. Counted row-wise, so a11, a12, a13, a21... etc
 Returns: the same, but inverted.
*/
void matrix_invert_3x3(float m[9], OUT float inverted_matrix[9]){
    float det = 1/(m[0]*m[4]*m[8]-m[0]*m[5]*m[7]-m[1]*m[3]*m[8]+m[1]*m[5]*m[6]+m[2]*m[3]*m[7]-m[2]*m[4]*m[6]);
    inverted_matrix[0] = (m[4]*m[8]-m[5]*m[7])*det;
    inverted_matrix[1] = (m[2]*m[7]-m[1]*m[8])*det;
    inverted_matrix[2] = (m[1]*m[5]-m[2]*m[4])*det;
    
    inverted_matrix[3] = (m[5]*m[6]-m[3]*m[8])*det;
    inverted_matrix[4] = (m[0]*m[8]-m[2]*m[6])*det;
    inverted_matrix[5] = (m[2]*m[3]-m[0]*m[5])*det;
    
    inverted_matrix[6] = (m[3]*m[7]-m[4]*m[6])*det;
    inverted_matrix[7] = (m[1]*m[6]-m[0]*m[7])*det;
    inverted_matrix[8] = (m[0]*m[4]-m[1]*m[3])*det;

}
// note: this needs to be run AFTER positioning
/**
 Get the equilibrium configuration of the attachment element.
 This uses the inverse Jacobian of the element to reconstruct the element
 at equilibrium. This is okay, because we don't need the absolute position
 of anything, only the relative positions.
 Params:
 * J_inv_0, the inverse jacobian of the tetrahedron.
 * ends_at_rod - if true, the attachment goes blob-to-rod, else it goes rod-to-blob
 * node_weighting - the positioning of the attachment node inside the rod, given as a fraction of the tetrahedron edges
 * tet_origin - origin of the tetrahedron (where the edge vectors originate)
 * edge_vecs - the edge vectors themselves.
 Returns:
 * equil_attachment_node - the equilibrium attachment element
*/
void equil_attachment_node_from_J(float J_inv_0[9], int face_node_indices[3], bool ends_at_rod, float node_weighting[3], float tet_origin[3], float edge_vecs[3][3], float rotation[3], OUT float equil_attachment_node[3]){
    float J[9];
    // reconstruct nodes from jacobian
    matrix_invert_3x3(J_inv_0, J);
    print_array("J", J, 9);
    float face_nodes[4][3];
    vec3d(n){face_nodes[0][n] = 0;}
    vec3d(n){face_nodes[1][n] = J[n];}
    vec3d(n){face_nodes[2][n] = J[n+3];}
    vec3d(n){face_nodes[3][n] = J[n+6];}
    
    print_array("face_node_1", face_nodes[0], 3);
    print_array("face_node_2", face_nodes[1], 3);
    print_array("face_node_3", face_nodes[2], 3);
    print_array("face_node_4", face_nodes[3], 3);
    
    // get attachment element as in get_attachment_node etc
    get_tri_norm(face_nodes[face_node_indices[0]], face_nodes[face_node_indices[1]], face_nodes[face_node_indices[2]], equil_attachment_node);
    normalize(equil_attachment_node, equil_attachment_node);
    
    float equil_attachment_node_pos_relative[3];
    get_attachment_node_pos(face_nodes[face_node_indices[0]], face_nodes[face_node_indices[1]], face_nodes[face_node_indices[2]], edge_vecs, node_weighting, tet_origin, equil_attachment_node_pos_relative);
    
    // set the direction of this element (accounting for the direction/'polarity' of the rod)
    bool points_out = points_out_of_tet(face_nodes[0], face_nodes[1], face_nodes[2], face_nodes[3], equil_attachment_node, equil_attachment_node_pos_relative);
    
    float polarity_multiplication_factor = -1;
    if (ends_at_rod){
        polarity_multiplication_factor = 1;
    }
    
    //loat dotprod = path_to_centroid[0]*attachment_node[0] + path_to_centroid[1]*attachment_node[1] + path_to_centroid[2]*attachment_node[2];
    if(dbg_print){std::cout << "   points_out" << points_out << "\n";}
    if (!points_out){
        if(dbg_print){std::cout << "   reversing attachment node.\n";}  //dbg
        vec3d(n){equil_attachment_node[n]*=polarity_multiplication_factor;}
    }
    
   float empty[3] = {0,0,0};
   if (array_equal(rotation, empty) == false){
       float euler_rm[9];
       get_rotation_matrix_from_euler(rotation, euler_rm);
       float rotated_attachment_node[3];
       apply_rotation_matrix(equil_attachment_node, euler_rm, rotated_attachment_node);
       vec3d(n){equil_attachment_node[n] = rotated_attachment_node[n];}
   }
    
}

/**
 Check to see whether an attachment element points out of or into its
 attached tetrahedron. If it points into the tetrahedron, it should be
 flipped around.
 Params:
 * node1, node2, node3, node4 - nodes in the tetrahedron
 * attachment_element - the attachment element
 * attachment_node - position of the attachment node
 Returns:
 * points_out, a boolean. True if it points out.
 
*/
bool points_out_of_tet(float node1[3], float node2[3], float node3[3], float node4[3], float attachment_element[3], float attachment_node[3]){
    float tet_centroid[3] = {0,0,0};
    vec3d(n){tet_centroid[n] += node1[n];}
    vec3d(n){tet_centroid[n] += node2[n];}
    vec3d(n){tet_centroid[n] += node3[n];}
    vec3d(n){tet_centroid[n] += node4[n];}
    vec3d(n){tet_centroid[n] /= 4.0;}
    float path_to_centroid[3];
    vec3d(n){path_to_centroid[n] = attachment_node[n] - tet_centroid[n];}
    float dotprod = 0;
    vec3d(n){dotprod += attachment_element[n] * path_to_centroid[n];}
    bool points_out = dotprod<0;
    return points_out;
}

/**
 equivalent to Rod::get_attachment_node_pos. Partial refactoring. Todo:
 replace the get_attachment_node_pos stuff inside rod_blob_interface
 with calls to this.
*/
void get_attachment_node_pos(float face_node_1[3], float face_node_2[3], float face_node_3[3], float edge_vecs[3][3], float node_weighting[3], float tet_origin[3], OUT float attachment_node_pos[3]){
    if (node_weighting[0] == -1 && node_weighting[1] == -1 && node_weighting[2] == -1){
        attachment_node_pos[0] = .33 * (face_node_1[0] + face_node_2[0] + face_node_3[0]);
        attachment_node_pos[1] = .33 * (face_node_1[1] + face_node_2[1] + face_node_3[1]);
        attachment_node_pos[2] = .33 * (face_node_1[2] + face_node_2[2] + face_node_3[2]);
    }
    else{ // otherwise use the weighting
        vec3d(n){attachment_node_pos[n] = tet_origin[n] + (edge_vecs[0][n]*node_weighting[0] + edge_vecs[1][n]*node_weighting[1] + edge_vecs[2][n]*node_weighting[2]);}
    }
}

// interface structure

/**
 Rod_blob_interface constructor. Requires pointers to fully-initialized rod and blob objects. Parameters:
 - set_ends_at_rod: if true, the connection goes from blob to rod. If false, it goes from rod to blob. This affects the index of the rod element which is connected.
 - set_to_index: index of the source element (via blob->get_element).
 - set_from_index: index of the destination element (via blob->get_element).
 - set_blob_node_ids: the indices of the three linear nodes on the blob which make up the face that the interface is orientated relative to.
 - rotation: an euler rotation matrix, which determines the rotation of the attachment element relative to the face it's attached to.
 - node_weighting: how much each of the first three edge vectors is weighted to determine the position of the attachment node inside the tetrahedron. You can also set them all to -1 to make the thing just go in the center of the element.
 - order: I think this is just here for reference, but it refers to the order the connection was resolved in (e.g. in a system of multiple connections, which one gets set up first)
 Anyway this is a boring constructor, all it does is set member variables, 
 allocate a bit of memory for the internal tetrahedron, set the
 initial internal state (that's the edge vectors and internal tetrahedron)
 and set the indices of the face nodes relative to the tetrahedron.
*/
Rod_blob_interface::Rod_blob_interface(Rod* set_connected_rod, Blob* set_connected_blob, bool set_ends_at_rod, int set_to_index, int set_from_index, int blob_node_ids[3], float rotation[3], float node_weighting[3], int order)
    {
        if(dbg_print){std::cout << "Creating rod.\n";}
        this->connected_rod = set_connected_rod;
        this->connected_blob = set_connected_blob;
        this->ends_at_rod = set_ends_at_rod;
        this->to_index = set_to_index;
        this->from_index = set_from_index;
        vec3d(n){this->face_nodes[n] = blob_node_ids[n];}
        this->order = order;
        vec3d(n){this->euler_angles[n] = rotation[n];}
        vec3d(n){this->node_weighting[n] = node_weighting[n];}
        
        if(dbg_print){std::cout << "  ends at rod: " << this->ends_at_rod << "\n";}
        print_array("  euler angles", this->euler_angles, 3);
        print_array("  node_weighting", this->node_weighting, 3);
        
        int element_no = this->get_element_id(blob_node_ids);
        if(dbg_print){std::cout << "element no: " << element_no << "\n";}
        
        if (ends_at_rod){
            //this->connected_tet = this->connected_blob->get_element(this->to_index);
            this->connected_tet = this->connected_blob->get_element(element_no);
        }
        else{
            this->connected_tet = this->connected_blob->get_element(element_no);
            //this->connected_tet = this->connected_blob->get_element(this->from_index);
        }
        
        if(dbg_print){std::cout << "  Element node indices: " << this->connected_tet->n[0]->index << ", " << this->connected_tet->n[1]->index << ", " << this->connected_tet->n[2]->index << ", " << this->connected_tet->n[3]->index << " | " << this->connected_tet->n[4]->index << ", " << this->connected_tet->n[5]->index << ", " << this->connected_tet->n[6]->index << ", " << this->connected_tet->n[7]->index << ", " << this->connected_tet->n[8]->index << ", " << this->connected_tet->n[9]->index << "\n";}
        
        //this->connected_face = this->connected_blob->absolutely_get_face(this->face_index*4); // second order faces are numbered oddly
        
        //vec3d(n){this->face_nodes[n] = this->connected_face->n[n]->index;}
        
        if(dbg_print){std::cout << "  Face nodes: " << this->face_nodes[0] << ", " << this->face_nodes[1] << ", " << this->face_nodes[2] << "\n";}
        if(dbg_print){std::cout << "  To index = " << this->to_index << "\n";}
        if(dbg_print){std::cout << "  From index = " << this->from_index << "\n";}
        
        //if (this->connected_face->e != this->connected_tet){
        //    std::cout << "To index = " << this->to_index << "\n";
        //    std::cout << "From index = " << this->from_index << "\n";
        //    rod_abort("The face at the rod_blob_interface does not belong to the tetrahedron of the rod_blob_interface.");
        //}
        
        for(int i =0; i<4; i++){
            this->deformed_tet_nodes[i] = new mesh_node();
        }

        this->update_internal_state(true, true);
        this->select_face_nodes(this->face_node_indices);
        
        //this->get_attachment_node(this->attachment_node_equil, this->attachment_node_pos_equil, true); // note: remove attachment_node_pos_equil, it is useless
        
        if (set_ends_at_rod){
            this->connected_rod->interface_at_start = true; //blob to rod
        }
        
        else{
            this->connected_rod->interface_at_end = true; //rod to blob
        }
        
        if(dbg_print){std::cout << "  Mesh node indices: " << this->connected_tet->n[0]->index << ", " << this->connected_tet->n[1]->index << ", " << this->connected_tet->n[2]->index << ", " << this->connected_tet->n[3]->index << "\n";}
        if(dbg_print){std::cout << "Rod created.\n";}
    }

void Rod_blob_interface::set_initial_values(){
        //matrix3 old_jacobian;
        //matrix3 old_jacobian_inverse;
        //scalar duuh;
        //this->connected_tet->calculate_jacobian(old_jacobian);
        //
        //if(dbg_print){std::cout << "J before : [" << old_jacobian[0][0] << ", " << old_jacobian[0][1] << ", " << old_jacobian[0][2] << ", " << old_jacobian[1][0] << ", " << old_jacobian[1][1] << ", " << old_jacobian[1][2] << ", " << old_jacobian[2][0] << ", " << old_jacobian[2][1] << ", " << old_jacobian[2][2] << "\n";}
        //
        //mat3_invert(old_jacobian, old_jacobian_inverse, &duuh);
        //for (int i=0; i<9; i++){
        //    int j = i/3;
        //    int k = i%3;
        //    this->J_inv_0[i] = (float)old_jacobian_inverse[j][k];
        //}
        
        for (int i=0; i<9; i++){
            int j = i/3;
            int k = i%3;
            this->J_inv_0[i] = (float)this->connected_tet->J_inv_0[j][k];
        }
        
        print_array("  J_inv_0", this->J_inv_0, 9);
        
        equil_attachment_node_from_J(this->J_inv_0, this->face_node_indices, this->ends_at_rod, this->node_weighting, this->tet_origin, this->edge_vecs, this->euler_angles, this->attachment_node_equil); // note: remove attachment_node_pos_equil, it is useless
        this->get_attachment_material_axis(this->attachment_node_equil, this->attachment_m_equil);
        
        this->get_attachment_node(this->attachment_node, this->attachment_node_pos, false);
        this->get_attachment_material_axis(this->attachment_node, this->attachment_m);
}

/**
 Update the internal state of the rod-blob interface. The update_edge_vecs parameter will reconstruct the edges, which are need for certain mathematical operations, including updating the position of the attachment node. update_tet updates the cached tetrahedron object stored by the interface.
 
 If you've deformed the tetrahedron somehow, you should leave this alone until you've computed the gradient deformation matrix, as you need 'before' and 'after' tetrahedra to make that work. If you're doing some numerical integraion, you should use the rod_blob_interface tetrahedron as the one you deform, because you can just update_tet as soon as you're done with it.
*/
void Rod_blob_interface::update_internal_state(bool update_edge_vecs, bool update_tet){
    if (update_tet){this->set_tet(this->connected_tet);}
    if (update_edge_vecs){this->set_edge_vecs();}
    //std::cout << "Updated internal state. \n"; \\dbg
    //print_array("  Edge vecs node 0", this->edge_vecs[0], 3); \\dbg
    //print_array("  Edge vecs node 1", this->edge_vecs[1], 3); \\dbg
    //print_array("  Edge vecs node 2", this->edge_vecs[2], 3); \\dbg
    //print_array("  deformed tetrahedron node 0", this->deformed_tet_nodes[0]->pos.data, 3); \\dbg
    //print_array("  deformed tetrahedron node 1", this->deformed_tet_nodes[1]->pos.data, 3); \\dbg
    //print_array("  deformed tetrahedron node 2", this->deformed_tet_nodes[2]->pos.data, 3); \\dbg
    //print_array("  deformed tetrahedron node 3", this->deformed_tet_nodes[3]->pos.data, 3); \\dbg
}

void Rod_blob_interface::update_J_0(){
    // only run this in initialisation, after box positioning, before restart loading
    this->update_internal_state(true, true);
    float J[9];
    get_jacobian(this->deformed_tet_nodes, J);
    float J_inv[9];
    matrix_invert_3x3(J, J_inv);
    for(int i=0; i<9; i++){this->J_inv_0[i] = J_inv[i];}
}

/**
 Recompute te edge vectors for a tetrahedron. These are not automatically updated when the nodes are moved. If you do something to the attachment tetrahedron, consider running Rod_blob_interface->update_internal_state().
*/
void Rod_blob_interface::set_edge_vecs(){
    vec3d(v){this->tet_origin[v] = this->deformed_tet_nodes[0]->pos.data[v];}
    vec3d(v){this->edge_vecs[0][v] = this->deformed_tet_nodes[1]->pos[v] - this->tet_origin[v];}
    vec3d(v){this->edge_vecs[1][v] = this->deformed_tet_nodes[2]->pos[v] - this->tet_origin[v];}
    vec3d(v){this->edge_vecs[2][v] = this->deformed_tet_nodes[3]->pos[v] - this->tet_origin[v];}
    
    //vec3d(v){this->tet_origin_equil[v] = this->deformed_tet_nodes[0]->pos_0.data[v];}
    //vec3d(v){this->edge_vecs_equil[0][v] = this->deformed_tet_nodes[1]->pos_0[v] - this->tet_origin_equil[v];}
    //vec3d(v){this->edge_vecs_equil[1][v] = this->deformed_tet_nodes[2]->pos_0[v] - this->tet_origin_equil[v];}
    //vec3d(v){this->edge_vecs_equil[2][v] = this->deformed_tet_nodes[3]->pos_0[v] - this->tet_origin_equil[v];}
}

/**
 Get the position and direction vectors for the current attachment node\element.
 Note: this function is great for initialisation, but once the connection
 is initialised, there's no need to run this. Instead, use 'reorientate connection'.
 The reason is that this function sets the material axis arbitrarily, whereas
 reorientate_connection will update the material axis based on the rotation
 of the tetrahedorn.
*/
void Rod_blob_interface::get_attachment_node(OUT float attachment_node[3], float attachment_node_pos[3], bool equil){ //equil
    
    // std::cout << "Getting attachment node\n"; //dbg
    
    //Note: do not look inside the function 'select_face_nodes'.
    //int face_node_indices[3];
    //this->select_face_nodes(face_node_indices);
    // std::cout << "  Face node indices: " << this->face_node_indices[0] << ", " << this->face_node_indices[1] << ", " << this->face_node_indices[2] << "\n"; //dbg
    
    
    //grab face nodes
    //vector3 face_node_1_data;
    //vector3 face_node_2_data;
    //vector3 face_node_3_data;
    //float face_node_1[3];
    //float face_node_2[3];
    //float face_node_3[3];
    //for (int i=0; i<3; i++){
    //    face_node_1[i] = this->connected_face->n[0]->pos[i];
    //    face_node_2[i] = this->connected_face->n[1]->pos[i];
    //    face_node_3[i] = this->connected_face->n[2]->pos[i];
    //    this->connected_blob->get_node(this->face_nodes[0], face_node_1_data.data);
    //    this->connected_blob->get_node(this->face_nodes[1], face_node_2_data.data);
    //    this->connected_blob->get_node(this->face_nodes[2], face_node_3_data.data);
    //}
    //vec3d(n){face_node_1[n] = face_node_1_data.data[n];}
    //vec3d(n){face_node_2[n] = face_node_2_data.data[n];}
    //vec3d(n){face_node_3[n] = face_node_3_data.data[n];}
    
    
    //rod::print_array("face node 1", face_node_1, 3);
    //rod::print_array("face node 2", face_node_2, 3);
    //rod::print_array("face node 3", face_node_3, 3);
    
    float face_node_1[3];
    float face_node_2[3];
    float face_node_3[3];
    float non_face_node[3];
    
    //if (equil){
    //    vec3d(n){face_node_1[n] = this->deformed_tet_nodes[this->face_node_indices[0]]->pos_0.data[n];}
    //    vec3d(n){face_node_2[n] = this->deformed_tet_nodes[this->face_node_indices[1]]->pos_0.data[n];}
    //    vec3d(n){face_node_3[n] = this->deformed_tet_nodes[this->face_node_indices[2]]->pos_0.data[n];}
    //}
    //else{
    vec3d(n){face_node_1[n] = this->deformed_tet_nodes[this->face_node_indices[0]]->pos.data[n];}
    vec3d(n){face_node_2[n] = this->deformed_tet_nodes[this->face_node_indices[1]]->pos.data[n];}
    vec3d(n){face_node_3[n] = this->deformed_tet_nodes[this->face_node_indices[2]]->pos.data[n];}
    //}
    
    for( int i=0; i<4; i++){
        bool i_in_face_nodes = false;
        for (int j = 0; j<3; j++){ if (face_node_indices[j] == i){i_in_face_nodes = true;} }
        if (!i_in_face_nodes){ vec3d(n){ non_face_node[n] = this->deformed_tet_nodes[i]->pos.data[n]; } }
    }
    
    if(dbg_print){std::cout << "face node indices: " << this->face_node_indices[0] << ", " << this->face_node_indices[1] << ", " << this->face_node_indices[2] << "\n";}
    
    print_array("face_node_1", face_node_1, 3);
    print_array("face_node_2", face_node_2, 3);
    print_array("face_node_3", face_node_3, 3);
    
    //get normal
    if (equil){
        equil_attachment_node_from_J(this->J_inv_0, this->face_node_indices, this->ends_at_rod, this->node_weighting, this->tet_origin, this->edge_vecs, this->euler_angles, attachment_node);
    }
    else{
        get_tri_norm(face_node_1, face_node_2, face_node_3, attachment_node);
    }
    
    print_array("attachment node", attachment_node, 3);
    
    normalize(attachment_node, attachment_node);
    
    print_array("normalized attachment node", attachment_node, 3);
    
    // note: if there is a rotation of the node needed, put it here!
    
    // print_array(" tet origin", this->tet_origin, 3); //dbg
    
    //get position of node in face
    //float minus1array[3] = {-1,-1,-1}; //set everything to -1 just to center the node
    this->update_internal_state(true, false);
    if (this->node_weighting[0] == -1 && this->node_weighting[1] == -1 && this->node_weighting[2] == -1){
        attachment_node_pos[0] = .33 * (face_node_1[0] + face_node_2[0] + face_node_3[0]);
        attachment_node_pos[1] = .33 * (face_node_1[1] + face_node_2[1] + face_node_3[1]);
        attachment_node_pos[2] = .33 * (face_node_1[2] + face_node_2[2] + face_node_3[2]);
        //attachment_node_pos[1] = .25 * (this->deformed_tet_nodes[0]->pos.y + this->deformed_tet_nodes[1]->pos.y + this->deformed_tet_nodes[2]->pos.y + this->deformed_tet_nodes[3]->pos.y);
        //attachment_node_pos[2] = .25 * (this->deformed_tet_nodes[0]->pos.z + this->deformed_tet_nodes[1]->pos.z + this->deformed_tet_nodes[2]->pos.z + this->deformed_tet_nodes[3]->pos.z);
    }
    else{ // otherwise use the weighting
        
        //if(equil){
        //vec3d(n){attachment_node_pos[n] = this->tet_origin_equil[n] + (this->edge_vecs_equil[0][n]*this->node_weighting[0] + this->edge_vecs_equil[1][n]*this->node_weighting[1] + this->edge_vecs_equil[2][n]*this->node_weighting[2]);}
        //}
        //else{
        vec3d(n){attachment_node_pos[n] = this->tet_origin[n] + (this->edge_vecs[0][n]*this->node_weighting[0] + this->edge_vecs[1][n]*this->node_weighting[1] + this->edge_vecs[2][n]*this->node_weighting[2]);}
       // }
    }

    float end_node[3];
    int index;
    if (this->ends_at_rod){
        index = 0;
    }
    else{
        index = this->connected_rod->num_nodes-1;
    }
    
    float end_node_pos[3];
    vec3d(n){end_node_pos[n] = this->connected_rod->current_r[(index*3)+n];}
    //normalize(end_node, end_node);
    
    print_array("   attachment_node_pos", attachment_node_pos, 3);
    print_array("   attachment_node", attachment_node, 3);
    //print_array("   centroid", centroid, 3);
    //print_array("   path to centroid", path_to_centroid, 3);

    // if we're going rod-to-blob, the rod is going into the blob, so the attachment element must point
    // toward the centroid. if we're going blob-to-rod, it must point out of the blob, so away from the centroid
    float polarity_multiplication_factor = -1;
    if (this->ends_at_rod){
        polarity_multiplication_factor = 1;
    }
    
    //float dotprod = path_to_centroid[0]*attachment_node[0] + path_to_centroid[1]*attachment_node[1] + path_to_centroid[2]*attachment_node[2];
    
    bool points_out = points_out_of_tet(face_node_1, face_node_2, face_node_3, non_face_node, attachment_node, attachment_node_pos);
    
    if(dbg_print){std::cout << "points out of tet: " << points_out << "\n";}
    
    //if(dbg_print){std::cout << "   dotprod" << dotprod << "\n";}
    //if (dotprod < 0){
    if (!points_out){
        if(dbg_print){std::cout << "   reversing attachment node.\n";}  //dbg
        vec3d(n){attachment_node[n]*=polarity_multiplication_factor;}
    }
    
    print_array("   attachment_node (reverse'd)", attachment_node, 3);
    vec3d(n){end_node[n] = this->connected_rod->current_r[(index*3)+3+n] - this->connected_rod->current_r[(index*3)+n];}
    print_array("   end node (unmodified)", end_node, 3);
    print_array("   end node pos (unmodified)", end_node_pos, 3);
        
    // print_array("  end_node", end_node, 3); //dbg

    //make sure normal is facing the right way
    //note: blob centroid, not the tet centroid
    //float way_to_centroid[3];
    //vector3 blob_centroid[3];
    //this->connected_blob->get_centroid(blob_centroid);
    //vec3d(n){way_to_centroid[n] = blob_centroid->data[n] - attachment_node_pos[n];}
    //vec3d(n){way_to_centroid[n] = this->connected_tet->centroid[n] - attachment_node_pos[n];}
    //normalize(way_to_centroid, way_to_centroid);
    //float dotprod = way_to_centroid[0]*attachment_node[0] + way_to_centroid[1]*attachment_node[1] + way_to_centroid[2]*attachment_node[2];
    // print_array("  attachment node pos", attachment_node_pos, 3); //dbg
    // print_array("  attachment node", attachment_node, 3); //dbg
    //print_array(" blob centroid", blob_centroid->data, 3);
    //print_array(" way_to_centroid", way_to_centroid, 3);
    // std::cout << "  dotprod: " << dotprod << "\n"; //dbg
    //if (dotprod < 0){
    //    std::cout << " reversing attachment node.\n";
    //    vec3d(n){attachment_node[n]*=-1;}
   // }
   
   float empty[3] = {0,0,0};
   if (array_equal(this->euler_angles, empty) == false){
       std::cout << "attachment node before rotation: " << attachment_node[0] << ", " << attachment_node[1] << ", " << attachment_node[2] << "\n";
       float euler_rm[9];
       get_rotation_matrix_from_euler(this->euler_angles, euler_rm);
       float rotated_attachment_node[3];
       apply_rotation_matrix(attachment_node, euler_rm, rotated_attachment_node);
       vec3d(n){attachment_node[n] = rotated_attachment_node[n];}
       std::cout << "attachment node after rotation: " << attachment_node[0] << ", " << attachment_node[1] << ", " << attachment_node[2] << "\n";
   }
   
   std::cout << "Debug: acquisition of interface.\n";
   std::cout << "Euler angles: " << this->euler_angles[0] << ", " << this->euler_angles[1] << ", " << this->euler_angles[2] << "\n";
   std::cout << "Points out: " << points_out << "\n";
   std::cout << "This interface order: " << this->order << "\n";
   
}
/**
 Get the position vector for the attachment node. Absolute, not relative
 to the tetrahedron or anything. Uses the weighting that the object
 was initialised with. 
*/
void Rod_blob_interface::get_attachment_node_pos(float attachment_node_pos[3], bool equil){

    this->update_internal_state(true, false);

    float face_node_1[3];
    float face_node_2[3];
    float face_node_3[3];
    
    if (equil){
        vec3d(n){face_node_1[n] = this->deformed_tet_nodes[this->face_node_indices[0]]->pos_0.data[n];}
        vec3d(n){face_node_2[n] = this->deformed_tet_nodes[this->face_node_indices[1]]->pos_0.data[n];}
        vec3d(n){face_node_3[n] = this->deformed_tet_nodes[this->face_node_indices[2]]->pos_0.data[n];}
    }
    else{
        vec3d(n){face_node_1[n] = this->deformed_tet_nodes[this->face_node_indices[0]]->pos.data[n];}
        vec3d(n){face_node_2[n] = this->deformed_tet_nodes[this->face_node_indices[1]]->pos.data[n];}
        vec3d(n){face_node_3[n] = this->deformed_tet_nodes[this->face_node_indices[2]]->pos.data[n];}
    }    

    
    if (this->node_weighting[0] == -1 && this->node_weighting[1] == -1 && this->node_weighting[2] == -1){
        attachment_node_pos[0] = 1./3. * (face_node_1[0] + face_node_2[0] + face_node_3[0]);
        attachment_node_pos[1] = 1./3. * (face_node_1[1] + face_node_2[1] + face_node_3[1]);
        attachment_node_pos[2] = 1./3. * (face_node_1[2] + face_node_2[2] + face_node_3[2]);
    }
    else{ // otherwise use the weighting
        
        //if(equil){
        //vec3d(n){attachment_node_pos[n] = this->tet_origin_equil[n] + (this->edge_vecs_equil[0][n]*this->node_weighting[0] + this->edge_vecs_equil[1][n]*this->node_weighting[1] + this->edge_vecs_equil[2][n]*this->node_weighting[2]);}
        //}
        //else{
        vec3d(n){attachment_node_pos[n] = this->tet_origin[n] + (this->edge_vecs[0][n]*this->node_weighting[0] + this->edge_vecs[1][n]*this->node_weighting[1] + this->edge_vecs[2][n]*this->node_weighting[2]);}
        //}
    }

}

// need to call this AFTER getting the attachment node
/**
 Once the attachment node has been obtained, this sets the initial direction of the attachment material axis.
 Note that the initial direction is arbitray - by default it is set such that there is no energy between the attachment axis and the first rod material axis.
 Only run this function for initialisation! When the simulation is running, you need to update the orientation of this thing using this->reorientate_connection.
*/
void Rod_blob_interface::get_attachment_material_axis(float attachment_node[3], OUT float attachment_material_axis[3]){
    float nearest_material_axis[3];
    float nearest_element[3];
    
    if (this->ends_at_rod){
        vec3d(n){ nearest_material_axis[n] = this->connected_rod->equil_m[this->connected_rod->length - 3 + n]; }
        vec3d(n){ nearest_element[n] = (this->connected_rod->equil_r[this->connected_rod->length - 3 + n]) - (this->connected_rod->equil_r[this->connected_rod->length - 6 + n]); }
    }
    
    else{
        vec3d(n){ nearest_material_axis[n] = this->connected_rod->equil_m[n]; }
        vec3d(n){ nearest_element[n] = this->connected_rod->equil_r[3 + n] - this->connected_rod->equil_r[n]; }
    }
    // the parallel transport algorithm is not stable for very specific configurations, including one
    // created by a particular (but common) set of rod-blob attachment parameters
    normalize(nearest_material_axis, nearest_material_axis);
    normalize(nearest_element, nearest_element);
    if (nearest_element[0] == -attachment_node[0] && nearest_element[1] == -attachment_node[1] && nearest_element[2] == -attachment_node[2]){
        std::cout << "A parallel transport error was corrected in the initial rod-blob connection configuration.\n";
        vec3d(n){attachment_material_axis[n] = nearest_material_axis[n];}
        return;
    }
    
    parallel_transport(nearest_material_axis, attachment_material_axis, nearest_element, attachment_node);
    
}

/**
 Overwrite the internal tetrahedron nodes with the nodes from a tetra_element_linear object.
 Only run this after the connection has been reorientated
*/
void Rod_blob_interface::set_tet(tetra_element_linear *tet){
    vec3d(n){this->deformed_tet_nodes[0]->pos[n] = tet->n[0]->pos[n];}
    vec3d(n){this->deformed_tet_nodes[1]->pos[n] = tet->n[1]->pos[n];}
    vec3d(n){this->deformed_tet_nodes[2]->pos[n] = tet->n[2]->pos[n];}
    vec3d(n){this->deformed_tet_nodes[3]->pos[n] = tet->n[3]->pos[n];}
    bool garbage_tet = false;
    for (int i = 0; i < 4; i++){
        for (int j = 0; j<3; j++){
            if (this->deformed_tet_nodes[i]->pos[j] > 1e20 || this->deformed_tet_nodes[i]->pos[j] < -1e20){
                garbage_tet = true;
            }
        }
    }
    
    if (garbage_tet){
        std::cout << "Warning: tetrahdron is full of garbage\n";
    }
}

/**
 For a given attachment element and material axis, get a new material axis and attachment element based on the rotation matrix found from the QR decomposition of the gradient deformation matrix. Run this every timestep!
*/
void Rod_blob_interface::reorientate_connection(float attachment_element_orig[3], float attachment_material_axis_orig[3], OUT float new_attachment_element[3], float new_attachment_material_axis[3]){
    float gradient_deformation[9] = {0,0,0,0,0,0,0,0,0};
    float R[9];
    float Q[9];
    get_gradient_deformation(this->J_inv_0, this->deformed_tet_nodes, gradient_deformation);
    QR_decompose_gram_schmidt(gradient_deformation, Q, R);
    apply_rotation_matrix(attachment_element_orig, Q, new_attachment_element); // yes, really
    apply_rotation_matrix(attachment_material_axis_orig, Q, new_attachment_material_axis);
    print_array("  R", R, 9); //dbg
    print_array("  Q", Q, 9); //dbg
}

/**
 Initialisation function. This will set the position of the rod relative
 to the blob. For some reason, it does it in this really arcane way, in
 which it gets a rotation matrix to the attachment node, applies it
 to each node and then repositions each element at the end of the previous
 element. One parameter: use_equil, if true it'll do it for the
 equilibrium rod instead of the current rod. This member function is
 called from world.cpp during intialisation.
*/
void Rod_blob_interface::position_rod_from_blob(bool use_equil){ // default false
    
    // create a copy of current_r
    vector<float> current_r_rotated(this->connected_rod->length);
    for (int i=0; i < this->connected_rod->length; i++){
        if (use_equil == false){ current_r_rotated[i] = this->connected_rod->current_r[i]; }
        else{ current_r_rotated[i] = this->connected_rod->equil_r[i]; }
    } 
    
    // get the 1st element and work out the rotation matrix that lines it up with the attachment node
    float p_0[3];
    float p_0_rotated[3];
    float attachment_node[3];
    float attachment_node_pos[3];
    float rm[9];
    this->connected_rod->get_p(0, p_0, use_equil);
    
    this->get_attachment_node(OUT attachment_node, attachment_node_pos, false);
    
    if(dbg_print){std::cout << "POSITIONING ROD FROM BLOB!\n";} //dbg
    print_array("-  attachment_node", attachment_node, 3); //dbg
    print_array("-  attachment_node_pos", attachment_node_pos, 3); //dbg
    
    float p_0_scale = sqrt(p_0[0]*p_0[0] + p_0[1]*p_0[1] + p_0[2]*p_0[2]);
    vec3d(n){p_0[n] /= p_0_scale;}
    
    get_rotation_matrix(p_0, attachment_node, rm);
    apply_rotation_matrix(p_0, rm, p_0_rotated);
    
    float *m_src;
    
    float m_0[3];
    if (use_equil){
        m_src = this->connected_rod->equil_m;
    }
    else{
        m_src = this->connected_rod->current_m;
    }
    
    // same as above for material axis
    vec3d(n){m_0[n] = m_src[n];}
    normalize(m_0, m_0);
    float m_0_rotated[3];
    apply_rotation_matrix(m_0, rm, m_0_rotated);
    vec3d(n){m_src[n] = m_0_rotated[n];}
    
    vec3d(n){p_0_rotated[n] *= p_0_scale;}
    
    // work out translation needed to line up rod with interface element
    float translation[3];
    vec3d(n){translation[n] = attachment_node_pos[n] - current_r_rotated[n];}
    vec3d(n){current_r_rotated[n] += translation[n];}
    vec3d(n){current_r_rotated[n+3] += (translation[n] + p_0[n]);}
    
    print_array("-  translation", translation, 3); //dbg
    
    print_array("p_0", p_0, 3); //dbg
    if(dbg_print){std::cout << "p_0_scale = " << p_0_scale << "\n";} //dbg
    print_array("p_0 rotated", p_0_rotated, 3); //dbg
    print_array("attachment_node_pos", attachment_node_pos, 3); //dbg
    print_array("rm", rm, 9); //dbg

    float p[3];
    float p_rotated[3];    
    float p_scale;
    
    float m[3];
    float m_rotated[3];
    
    vec3d(n){current_r_rotated[n] = attachment_node_pos[n];}
    vec3d(n){current_r_rotated[n+3] = attachment_node_pos[n] + p_0_rotated[n];}
    
    // build the wretched thing one element at a time
    for (int i=1; i<this->connected_rod->num_nodes - 1; i++){
        //r
        this->connected_rod->get_p(i, p, use_equil);
        p_scale = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        vec3d(n){p[n] /= p_scale;}
        apply_rotation_matrix(p, rm, p_rotated);
        vec3d(n){p_rotated[n] *= p_scale;}
        print_array("p_rotated", p_rotated, 3);  //dbg
        vec3d(n){current_r_rotated[((i+1)*3)+n] = current_r_rotated[(i*3)+n] + p_rotated[n];}
        
        //m
        vec3d(n){m[n] = m_src[n+(i*3)];}
        normalize(m, m);
        apply_rotation_matrix(m, rm, m_rotated);
        vec3d(n){m_src[n+(i*3)] = m_rotated[n];}
    }
    
    for (int i=0; i < this->connected_rod->length; i++){
        if (use_equil == false){ this->connected_rod->current_r[i] = current_r_rotated[i]; }
        else{ this->connected_rod->equil_r[i] = current_r_rotated[i]; }
    } 
    
}

/**
 This one is fairly normal, in contrast to rod_form_blob. It gets the
 rotation matrix that aligns the attachment element with the rod element
 at the end of the rod, then applies that to the blob. Like
 position_rod_from_blob, it's called during initialisation from world.cpp.
*/
void Rod_blob_interface::position_blob_from_rod(){
    float p_end[3];
    float p_end_normalized[3];
    float attachment_node[3];
    float attachment_node_pos[3];
    float rm[9];
    float rod_end[3];
    float translate[3];
    this->update_internal_state(true, true);
    
    // get end node
    this->connected_rod->get_p(connected_rod->num_nodes-2, p_end, true); // -2! indexed from 0! duh
    get_attachment_node(attachment_node, attachment_node_pos, false);
    normalize(p_end, p_end_normalized);
    
    // get rotation matrix that aligns end node with attachment node, apply it to blob
    get_rotation_matrix(attachment_node, p_end_normalized, rm);
    this->connected_blob->rotate(rm[0], rm[1], rm[2], rm[3], rm[4], rm[5], rm[6], rm[7], rm[8], true );
    
    // the same as above but now for position instead of rotation
    this->update_internal_state(true, true);
    float attachment_node_pos_rotated[3];
    this->get_attachment_node_pos(attachment_node_pos_rotated, false);
    vec3d(n){rod_end[n] = this->connected_rod->current_r[this->connected_rod->length-3+n];}
    vec3d(n){translate[n] = rod_end[n] - attachment_node_pos_rotated[n];}
    this->connected_blob->move(translate[0], translate[1], translate[2]);
    this->update_internal_state(true, true);
    float attachment_node_pos_after[3];
    get_attachment_node_pos(attachment_node_pos_after, false);
    
    std::cout << "\n";
    float new_rod_end[3]; vec3d(n){new_rod_end[n] = this->connected_rod->current_r[this->connected_rod->length-3+n];}
    print_array("rod_end_before (should be some random ass thing)", rod_end, 3);
    print_array("attachment_node_pos_rotated", attachment_node_pos_rotated, 3);
    print_array("translate (rod_end - attachment_node_pos_rotated", translate, 3);
    print_array("attachment_node_pos_after (attachment_node_pos_rotated + translate)", attachment_node_pos_after, 3);
    std::cout << "\n";
    
    float new_attachment_node_pos[3]; 
    get_attachment_node_pos(new_attachment_node_pos, false);
    
    float new_attachment_node[3];
    get_attachment_node(new_attachment_node, attachment_node_pos, false);
    
    vec3d(n){this->attachment_node_equil[n] = new_attachment_node[n];}

}

/**
 Given the indices of the face nodes relative to the entire blob,
 this will give you the indices of the face nodes relative just to one
 tetrahedron. I say 'given indices of face nodes' but that's really given
 as a parameter in the constructor.
 Note: do not read the body of this function, just move on
*/
void Rod_blob_interface::select_face_nodes(OUT int face_node_indices[3]){
    
    vector3 face_node_1_data;
    vector3 face_node_2_data;
    vector3 face_node_3_data;
    float face_nodes[3][3];
    for (int i=0; i<3; i++){
        this->connected_blob->get_node(this->face_nodes[0], face_node_1_data.data);
        this->connected_blob->get_node(this->face_nodes[1], face_node_2_data.data);
        this->connected_blob->get_node(this->face_nodes[2], face_node_3_data.data);
    }
    vec3d(n){face_nodes[0][n] = face_node_1_data.data[n];}
    vec3d(n){face_nodes[1][n] = face_node_2_data.data[n];}
    vec3d(n){face_nodes[2][n] = face_node_3_data.data[n];}
    
    float all_nodes[4][3];
    for (int i=0; i<4; i++){
        vec3d(n){all_nodes[i][n] = (float)this->deformed_tet_nodes[i]->pos.data[n];}
    }
    
    // check for equality
    for (int i=0; i<4; i++){
        for (int j=0; j<3; j++){
            if(dbg_print){std::cout << "testing i=" << i << " and j=" << j << "\n";}
            print_array("First thing", all_nodes[i], 3);
            print_array("Second thing", face_nodes[j], 3);
            if (array_equal(all_nodes[i], face_nodes[j])){face_node_indices[j] = i;}
        }
    }
    
}

/**
 Similar to the above function, but it instead gets the ID of an element
 given 3 nodes inside it.
*/
int Rod_blob_interface::get_element_id(int nodes[3]){
    
    int element_id = -1;
    
    vector3 face_node_1_data;
    vector3 face_node_2_data;
    vector3 face_node_3_data;
    float face_nodes[3][3];
    for (int i=0; i<3; i++){
        this->connected_blob->get_node(this->face_nodes[0], face_node_1_data.data);
        this->connected_blob->get_node(this->face_nodes[1], face_node_2_data.data);
        this->connected_blob->get_node(this->face_nodes[2], face_node_3_data.data);
    }
    vec3d(n){face_nodes[0][n] = face_node_1_data.data[n];}
    vec3d(n){face_nodes[1][n] = face_node_2_data.data[n];}
    vec3d(n){face_nodes[2][n] = face_node_3_data.data[n];}
    
    // Note: Unfortunate as it is, this cannot be avoided.
    
    int num_elements = this->connected_blob->get_num_elements();
    for (int i=0; i<num_elements; i++){
        
        tetra_element_linear* curr_element;
        curr_element = this->connected_blob->get_element(i);
        float elem_nodes[4][3];
        
        for (int j=0; j<4; j++){
            vec3d(n){elem_nodes[j][n] = (float)curr_element->n[j]->pos.data[n];}
        }

        if( array_contains(elem_nodes, face_nodes)){
            element_id = i;
        }
    }
    
    return element_id;
}

/**
 This is the rod-blob interface equivalent of get_perturbation_energy. It
 does one-half of the dynamics necessary for the connection to work, namely
 it works out the forces being transferred from the rod to the nodes of the
 connected tetrahedron.
 Params:
  - node_index, the index (relative to the element) of the node that we're getting the energies for
  - attachment_node_equil - equilibrium position of the attachment node (when I say 'equilibrium' I mean 'with the tetrahedron in its current state, unaltered')
  - attachment_material_axis - same, but material axis
  - displacement - how much the node is being moved. I would suggest this->connected_rod->perturbation_amount
  - energy - this is the output listing the energy associated with that perturbation in the following axes: [+x +y +z -x -y -z].
*/
void Rod_blob_interface::get_node_energy(int node_index, float attachment_node_equil[3], float attachment_material_axis_equil[3], float attachment_node[3], float attachment_material_axis[3], float displacement, float energy[6]){ //int direction, float equil_energy, bool use_equil_energy){
    
    float attachment_node_pos[3];
    float attachment_n[3];
    float attachment_n_equil[3];
    float p[2][3];
    float m[2][3];
    float n[2][3];
    float p_equil[2][3];
    float m_equil[2][3];
    float n_equil[2][3];
    float k;
    float beta[2];
    float B[2][4];
    
    if(dbg_print){std::cout << "Getting node energy for node " << node_index << "\n";}
    
    //  set up indices
    //  adjacent index is the index of the thing next to the attachment node,
    //  double_adjacent_index is one node over, which we need in order to get the energies.
    //  these are different depending on which end of the rod the attachment is at.
    int adjacent_index;
    int double_adjacent_index;
    int triple_adjacent_index;
    
    int mat_adjacent_index;
    int mat_double_adjacent_index;
    int mat_triple_adjacent_index;
    
    // separate indices for material axes which are shifted backward one for ends_at_blob
    if (this->ends_at_rod){
        if(dbg_print){std::cout << " Ends at rod \n";}
        adjacent_index = 0;
        double_adjacent_index = 1;
        triple_adjacent_index = 2;
        mat_adjacent_index = 0;
        mat_double_adjacent_index = 1;
        mat_triple_adjacent_index = 2;

    }
    else{
        if(dbg_print){std::cout << " Ends at blob \n";}
        adjacent_index = this->connected_rod->num_nodes-1;
        double_adjacent_index = this->connected_rod->num_nodes-2;
        triple_adjacent_index = this->connected_rod->num_nodes-3;
        mat_adjacent_index = this->connected_rod->num_nodes-2;
        mat_double_adjacent_index = this->connected_rod->num_nodes-3;
        mat_triple_adjacent_index = this->connected_rod->num_nodes-4;
    }

    // get equil p, m
    // to get p, just connect the existing nodes together with a new element
    vec3d(n){p_equil[0][n] = this->connected_rod->equil_r[(double_adjacent_index*3)+n] - this->connected_rod->equil_r[(adjacent_index*3)+n];}
    vec3d(n){p_equil[1][n] = this->connected_rod->equil_r[(triple_adjacent_index*3)+n] - this->connected_rod->equil_r[(double_adjacent_index*3)+n];}

    // if ends_at_rod, the elements end up backwards, this flips them around (nothing needed for m)
    if (ends_at_rod == false){
        vec3d(n){p_equil[0][n] *= -1;}
        vec3d(n){p_equil[1][n] *= -1;}
        vec3d(n){p[0][n] *= -1;}
        vec3d(n){p[1][n] *= -1;}
    }
    
    // to get m, parallel transport the double adjacent node onto the new element we've created
    vec3d(n){m_equil[1][n] = this->connected_rod->equil_m[(mat_double_adjacent_index*3)+n];}
    vec3d(n){m_equil[0][n] = this->connected_rod->equil_m[(mat_adjacent_index*3)+n];}
    float p_0_equil_norm[3];
    float p_1_equil_norm[3];
    normalize(m_equil[1], m_equil[1]);
    normalize(m_equil[0], m_equil[0]);
    normalize(p_equil[0], p_0_equil_norm);
    normalize(p_equil[1], p_1_equil_norm);

    cross_product(attachment_material_axis_equil, attachment_node_equil, attachment_n_equil);
    cross_product(p_0_equil_norm, m_equil[0], n_equil[0]);
    cross_product(p_1_equil_norm, m_equil[1], n_equil[1]);

    // get constants
    // all rods are now expected to have constants for end nodes as these become the interface constants
    if (this->ends_at_rod){
        k = this->connected_rod->material_params[adjacent_index*3];
    }
    else{
        k = this->connected_rod->material_params[double_adjacent_index*3];
    }
    
    beta[0] = this->connected_rod->material_params[(adjacent_index*3)+1];
    beta[1] = this->connected_rod->material_params[(double_adjacent_index*3)+1];
    
    for (int i=0; i<4; i++){ B[0][i] = this->connected_rod->B_matrix[(adjacent_index*4)+i]; }
    for (int i=0; i<4; i++){ B[1][i] = this->connected_rod->B_matrix[(double_adjacent_index*4)+i]; }
    
    print_array(" undeformed tetrahedron node 0", this->deformed_tet_nodes[0]->pos.data, 3);
    print_array(" undeformed tetrahedron node 1", this->deformed_tet_nodes[1]->pos.data, 3);
    print_array(" undeformed tetrahedron node 2", this->deformed_tet_nodes[2]->pos.data, 3);
    print_array(" undeformed tetrahedron node 3", this->deformed_tet_nodes[3]->pos.data, 3);
    print_array(" attachment_node_equil", attachment_node_equil, 3);
    print_array(" attachment_n_equil", attachment_n_equil, 3);
    print_array(" attachment_material_axis_equil", attachment_material_axis_equil, 3);
    print_array(" p_equil[0]", p_equil[0], 3); //dbg
    print_array(" p_equil[1]", p_equil[1], 3); //dbg
    print_array(" m_equil[0]", m_equil[0], 3); //dbg
    print_array(" m_equil[1]", m_equil[1], 3); //dbg
    print_array(" n_equil[0]", n_equil[0], 3); //dbg
    print_array(" n_equil[1]", n_equil[1], 3); //dbg
    print_array(" B[0]", B[0], 4); //dbg
    print_array(" B[1]", B[1], 4); //dbg
    print_array(" beta", beta, 2); //dbg
    if(dbg_print){std::cout << " k:  " << k << "\n";} //dbg
    
    for(int i=0; i<6; i++){ //dimensions, + and -
        
        if(dbg_print){std::cout << "  i = " << i << "\n";}
        
        // set up iteration over node displacement in 3d
        int displacement_sign = -1;
        int displacement_index;
        if (i<3){ displacement_sign = 1; }
        if (i<3){
            displacement_index = i;
        }
        else{
            displacement_index = i-3;
        }
        
        if(dbg_print){std::cout << " ...with displacement " << displacement << " in axis " << displacement_index << " with sign " << displacement_sign << "\n";}
        
        // apply displacement to tetrahedron
        this->deformed_tet_nodes[node_index]->pos.data[displacement_index] += (displacement_sign*(displacement));
        
        // get actual attachment node
        this->reorientate_connection(attachment_node_equil, attachment_material_axis_equil, attachment_node, attachment_material_axis);
                
        float facevec[3];
        vec3d(n){facevec[n] = deformed_tet_nodes[this->face_node_indices[0]][n].pos.data -  deformed_tet_nodes[this->face_node_indices[1]][n].pos.data;}
        normalize(facevec, facevec);
        //float edge_attachment_dotprod = (facevec[0] * attachment_node[0]) + (facevec[1] * attachment_node[1]) + (facevec[2] * attachment_node[2]);
                
        this->get_attachment_node_pos(attachment_node_pos, false);
        
        this->position_rod_ends(attachment_node_pos);
        
        // get p,m
        vec3d(n){p[0][n] = this->connected_rod->current_r[(double_adjacent_index*3)+n] - this->connected_rod->current_r[(adjacent_index*3)+n];}
        vec3d(n){p[1][n] = this->connected_rod->current_r[(triple_adjacent_index*3)+n] - this->connected_rod->current_r[(double_adjacent_index*3)+n];}

        // if ends at rod the elements end up backwards, this flips them around
        if (ends_at_rod == false){
            vec3d(n){p[0][n] *= -1;}
            vec3d(n){p[1][n] *= -1;}
        }

        vec3d(n){m[1][n] = this->connected_rod->current_m[(mat_double_adjacent_index*3)+n];}
        vec3d(n){m[0][n] = this->connected_rod->current_m[(mat_adjacent_index*3)+n];}

        // normalize
        float p_0_norm[3];
        float p_1_norm[3];
        normalize(m[1], m[1]);
        normalize(m[0], m[0]);
        normalize(p[0], p_0_norm);
        normalize(p[1], p_1_norm);
        
        cross_product(attachment_material_axis, attachment_node, attachment_n);
        cross_product(p_0_norm, m[0], n[0]);
        cross_product(p_1_norm, m[1], n[1]);
        
        // energy lookup table!
        // energy layout = [+x, +y, +z, -x, -y, -z]
        energy[i] += get_stretch_energy(k, p[0], p_equil[0]);
        
        // attachment node is scaled to provide neutral weighting for bend energy
        rescale_attachment_node(attachment_node, p[0], attachment_node_equil, p_equil[0], attachment_node, attachment_node_equil);
        
        print_array(" deformed tetrahedron node 0", this->deformed_tet_nodes[0]->pos.data, 3);
        print_array(" deformed tetrahedron node 1", this->deformed_tet_nodes[1]->pos.data, 3);
        print_array(" deformed tetrahedron node 2", this->deformed_tet_nodes[2]->pos.data, 3);
        print_array(" deformed tetrahedron node 3", this->deformed_tet_nodes[3]->pos.data, 3);
        print_array(" attachment_node", attachment_node, 3);
        print_array(" attachment_node_pos", attachment_node_pos, 3);
        print_array(" attachment_material_axis", attachment_material_axis, 3);
        print_array(" attachment_n", attachment_n, 3);
        print_array(" p[0]", p[0], 3);
        print_array(" p[1]", p[1], 3);
        print_array(" n[0]", n[0], 3);
        print_array(" n[1]", n[1], 3);
        print_array(" m[0]", m[0], 3);
        print_array(" m[1]", m[1], 3);
        print_array(" p_equil[0]", p_equil[0], 3);
        print_array(" p_equil[1]", p_equil[1], 3);
        print_array(" n_equil[0]", n_equil[0], 3);
        print_array(" m_equil[0]", m_equil[0], 3);
        print_array(" n_equil[1]", n_equil[1], 3);
        print_array(" m_equil[1]", m_equil[1], 3);

        energy[i] += get_bend_energy_mutual_parallel_transport(attachment_node, p[0], attachment_node_equil, p_equil[0], attachment_n, attachment_material_axis, attachment_n_equil, attachment_material_axis_equil, n[0], m[0], n_equil[0], m_equil[0], B[0], B[0]);
        energy[i] += get_bend_energy_mutual_parallel_transport(p[0], p[1], p_equil[0], p_equil[1], n[0], m[0], n_equil[0], m_equil[0], n[1], m[1], n_equil[1], m_equil[1], B[1], B[1]);
        energy[i] += get_twist_energy(beta[0], m[0], attachment_material_axis, m_equil[0], attachment_material_axis_equil, attachment_node, p[0], attachment_node_equil, p_equil[0]);
        energy[i] += get_twist_energy(beta[1], m[1], m[0], m_equil[1], m_equil[0], p[0], p[1], p_equil[0], p_equil[1]);
        
        normalize(attachment_node, attachment_node);
        normalize(attachment_node_equil, attachment_node_equil);

        // un-apply displacement to tetrahedron
        this->deformed_tet_nodes[node_index]->pos.data[displacement_index] += (displacement_sign*displacement*-1);
        this->position_rod_ends(attachment_node_pos);
        
    }
    
    if(dbg_print){std::cout << "EXPLODING energy for node " << node_index << " = " << energy[0] << " " << energy[1] << " " << energy[2] << " " << energy[3] << " " << energy[4] << " " << energy[5] << "\n";}


}

/**
 Fix the position of the final rod node at the interface node position.
 Also, update the material axis associated with the end element, the
 same way as it's done in the rod dynamics (using update_m1_matrix).
 The energies about this end node are actually what transmit forces
 from the blob to the rod, so they're important!
*/
void Rod_blob_interface::position_rod_ends(float attachment_node_pos[3]){
        if (this->ends_at_rod){ //blob to rod
            float p_0[3]; float m_0[3];
            float p_0_prime[3]; float m_0_prime[3];
            vec3d(n){p_0[n] = this->connected_rod->current_r[3+n] - this->connected_rod->current_r[n];}
            vec3d(n){m_0[n] = this->connected_rod->current_m[n];}
            vec3d(n){p_0_prime[n] = this->connected_rod->current_r[3+n] - attachment_node_pos[n];}
            update_m1_matrix(m_0, p_0, p_0_prime, m_0_prime);
            vec3d(n){this->connected_rod->current_r[n] = attachment_node_pos[n];}
            vec3d(n){this->connected_rod->current_m[n] = m_0_prime[n];}
        }
        else{ //rod to blob
            float p_n[3]; float m_n[3];
            float p_n_prime[3]; float m_n_prime[3];
            vec3d(n){p_n[n] = this->connected_rod->current_r[this->connected_rod->length-6+n] - this->connected_rod->current_r[this->connected_rod->length-3+n];}
            vec3d(n){m_n[n] = this->connected_rod->current_m[this->connected_rod->length-6+n];}
            vec3d(n){p_n_prime[n] = attachment_node_pos[n] - this->connected_rod->current_r[this->connected_rod->length-6+n];}
            update_m1_matrix(m_n, p_n, p_n_prime, m_n_prime);
            vec3d(n){this->connected_rod->current_r[this->connected_rod->length-3+n] = attachment_node_pos[n];}
            vec3d(n){this->connected_rod->current_m[this->connected_rod->length-3+n] = m_n_prime[n];}
        }
}


/**
 This is the interface equivalent to the FFEA_rod 'do_timestep' function,
 and does largely the same things.
 1) Update the geometric state of the interface to account for how the blob has moved.
 2) Compute the energy gradient resulting from the interface at each tetrahedron node
 3) Compute the force on each tetrahedron node as if it were a rod node, and apply that force to the node
 No params or return values for this one, though you need to have a
 correctly initialized rod-blob interface, which means you have to have
 run set_initial_values. In regular FFEA dynamics, this happens after
 the dynamics of the rod/blob, though it's kind of arbitrary when to
 run it.
*/
void Rod_blob_interface::do_connection_timestep(){ // run this after regular blob\rod dynamics
    
    this->reorientate_connection(this->attachment_node, this->attachment_m, this->attachment_node, this->attachment_m);
    this->update_internal_state(true, true);
    this->get_attachment_node_pos(this->attachment_node_pos, false);
    this->position_rod_ends(attachment_node_pos);
    float dynamics_displacement = this->connected_rod->perturbation_amount;
    
    // For each tetrahedron node:
    for(int i=0; i<4; i++){
        float curr_node_energy[6] = {0,0,0,0,0,0};
        get_node_energy(i, this->attachment_node_equil, this->attachment_m_equil, this->attachment_node, this->attachment_m, dynamics_displacement*0.5, curr_node_energy)   ; 
        vector3 force;
        vec3d(n){force[n] = (curr_node_energy[n+3] - curr_node_energy[n])/dynamics_displacement;}
        print_array("Force added to node: ", force.data, 3);
        //std::cout << "Interface " << this->order << " node " << i << " force: [" << force.data[0] << ", " << force.data[1] << ", " << force.data[2] << "]\n"; //DEBUGGO
        this->connected_blob->add_force_to_node(force, this->connected_tet->n[i]->index);
    }
    
    this->update_internal_state(true, true);
    
}

}
