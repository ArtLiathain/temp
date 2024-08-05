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

#ifndef ROD_BLOB_INTERFACE
#define ROD_BLOB_INTERFACE

#include "rod_structure.h"

// this bit
//#include "SimulationParams.h"
#include "tetra_element_linear.h"
#include "Face.h"
#include "Blob.h"

namespace rod {

// functions go here

void get_tri_norm(float node0[3], float node1[3], float node2[3], OUT float tri_norm[3]);
void get_jacobian(mesh_node **tet_nodes, OUT float J[9]);
void float_3x3_invert(float m[9], OUT float m_inv[9]);
void get_gradient_deformation(float J_inv_0[9], mesh_node**nodes_curr, OUT float gradient_deformation_3x3[9]);
void QR_decompose_gram_schmidt(float matrix_3x3[9], OUT float Q[9], float R[9]);
void construct_euler_rotation_matrix(float a, float b, float g, float rotmat[9]);
void rotate_tet(float rotmat[9], mesh_node **nodes, OUT mesh_node **rotated_nodes);
void get_euler_angles(float rm[9], OUT float euler_angles[3]);
void get_rotation_matrix_from_euler(float euler_angles[3], OUT float rm[9]);
bool array_equal(float arr1[3], float arr2[3]);
bool array_contains(float large_arr[4][3], float small_arr[3][3]);
void mesh_node_null_check(mesh_node* node, std::string location);
void rescale_attachment_node(float attachment_node[3], float end_node[3], float attachment_node_equil[3], float end_node_equil[3], OUT float scaled_attachment_node[3], float scaled_attachment_node_equil[3]);
void matrix_invert_3x3(float m[9], OUT float inverted_matrix[9]);
void equil_attachment_node_from_J(float J_inv_0[9], int face_node_indices[3], bool ends_at_rod, float node_weighting[3], float tet_origin[3], float edge_vecs[3][3], float rotation[3], OUT float equil_attachment_node[3]);
bool points_out_of_tet(float node1[3], float node2[3], float node3[3], float node4[3], float attachment_element[3], float attachment_node[3]);
void get_attachment_node_pos(float face_node_1[3], float face_node_2[3], float face_node_3[3], float edge_vecs[3][3], float node_weighting[3], float tet_origin[3], OUT float face_node_pos[3]);

// objects go here

struct Rod_blob_interface
{
    // member variables
    Rod* connected_rod;
    Blob* connected_blob;
    Face* connected_face;
    tetra_element_linear* connected_tet;
    bool ends_at_rod = true; // if this is false, it goes rod->blob, otherwise it goes blob->rod
    int to_index;
    int from_index;
    int face_index;
    int order;
    int face_nodes[3];
    float node_weighting[3] = {0.333333333333333, 0.333333333333333, 0.333333333333333};
    float euler_angles[3] = {0, 0, 0};
    float tet_origin[3];
    //float tet_origin_equil[3];
    float edge_vecs[3][3];
    //float edge_vecs_equil[3][3];
    mesh_node* deformed_tet_nodes[4]; // used for calculating the jacobian
    int face_node_indices[3];
    float J_inv_0[9]; // equilibrium jacobian of the attachment node
    
    float attachment_node_equil[3];
    //float attachment_node_pos_equil[3]; // note: remove attachment_node_pos_equil, it is useless
    float attachment_m_equil[3];
    
    float attachment_node[3];
    float attachment_node_pos[3];
    float attachment_m[3];
    
    //methods
    Rod_blob_interface (Rod* set_connected_rod, Blob* set_connected_blob, bool set_ends_at_rod, int set_to_index, int set_from_index, int face_nodes[3], float rotation[3], float node_weighting[3], int order);
    void set_initial_values();
    void update_internal_state(bool update_edge_vecs, bool update_tet);
    void update_J_0();
    void set_edge_vecs();
    void get_attachment_node(OUT float attachment_node[3], float attachment_node_pos[3], bool equil);
    void get_attachment_node_pos(float attachment_node_pos[3], bool equil);
    void get_initial_material_axis(); // parallel transport nearest mataxis
    void get_attachment_material_axis(float attachment_node[3], OUT float attachment_material_axis[3]);
    void get_tet_rotation_matrix(float tet_points_before[12], float tet_points_after[12]); // calls to get_gradient_deformation and qr_decompose
    
    void make_tet(); // turn nodes into tet elements
    void set_tet(tetra_element_linear *tet);
    void reorientate_connection(float attachment_element_orig[3], float attachment_material_axis_orig[3], OUT float new_attachment_element[3], float new_attachment_material_axis[3]);
    void position_rod_from_blob(bool use_equil);
    void position_blob_from_rod();
    void get_face(OUT float face_node_1[3], float face_node_2[3], float face_node_3[3]);
    void select_face_nodes(OUT int face_node_indices[3]);
    int get_element_id(int nodes[3]);
    void get_node_energy(int node_index, float attachment_node_equil[3], float attachment_material_axis_equil[3], float attachment_node[3], float attachment_material_axis[3], float displacement, float energy[6]);
    //void get_rod_energy(float attachment_node_equil[3], float attachment_material_axis_equil[3], float displacement, float energy[2][6]);
    void position_rod_ends(float attachment_node_pos[3]);
    void do_connection_timestep();
    

    // note: maybe a wrapper function for doubles that converts them to floats? see if it works
    
    // note: use vector12 class for internal nodes or something similar
    
};

bool point_inside_tetrahedron(float point[3], float tet[12]); // work out what direction the attachment node should face (it should point inside the tet, i think?)

} //end namespace


#endif
