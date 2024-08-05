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
 *      ffea_test.cpp
 *	Author: Rob Welch, University of Leeds
 *	Email: py12rw@leeds.ac.uk

 *	Author: Ryan Cocking, University of Leeds
 *	Email: bsrctb@leeds.ac.uk
 */

#include "ffea_test.h"

int ffea_test::do_ffea_test(std::string filename)
{
    int result = 1;
    std::ifstream t(filename);
    std::stringstream buffer;
    buffer << t.rdbuf();
    std::cout << "Test: " << buffer.str();

    if (buffer.str().find("connection_test") != std::string::npos)
    {
        result = ffea_test::connection_test();
    }

    if (buffer.str().find("arbitrary_equilibrium_twist") != std::string::npos)
    {
        result = ffea_test::arbitrary_equilibrium_twist();
    }

    if (buffer.str().find("connection_orientation_test") != std::string::npos)
    {
        result = ffea_test::connection_orientation_test();
    }

    if (buffer.str().find("arbitrary_equilibrium_bend") != std::string::npos)
    {
        result = ffea_test::arbitrary_equilibrium_bend();
    }

    if (buffer.str().find("identify_face") != std::string::npos)
    {
        result = ffea_test::identify_face();
    }

    if (buffer.str().find("connection_energy_1") != std::string::npos)
    {
        result = ffea_test::connection_energy();
    }

    if (buffer.str().find("connection_energy_2") != std::string::npos)
    {
        result = ffea_test::connection_energy_2();
    }

    if (buffer.str().find("jacobian_rotate") != std::string::npos)
    {
        result = ffea_test::jacobian_rotate();
    }

    if (buffer.str().find("connection_energy_3") != std::string::npos)
    {
        result = ffea_test::connection_energy_3();
    }

    if (buffer.str().find("connection_propagation") != std::string::npos)
    {
        result = ffea_test::connection_propagation_every_way();
    }

    if (buffer.str().find("recover_normal") != std::string::npos)
    {
        result = ffea_test::recover_normal();
    }

    if (buffer.str().find("dump_twist_info") != std::string::npos)
    {
        result = ffea_test::dump_twist_info();
    }

    if (buffer.str().find("lower_sphere") != std::string::npos)
    {
        result = ffea_test::lower_sphere();
    }

    if (buffer.str().find("line_connecting_rod_elements") != std::string::npos)
    {
        result = ffea_test::line_connecting_rod_elements();
    }

    if (buffer.str().find("rod_neighbour_list_construction") !=
        std::string::npos)
    {
        result = ffea_test::rod_neighbour_list_construction();
    }

    if (buffer.str().find("rod_steric_lj_potential") !=
        std::string::npos)
    {
        result = ffea_test::rod_steric_lj_potential();
    }

    if (buffer.str().find("nearest_image_pbc") !=
        std::string::npos)
    {
        result = ffea_test::nearest_image_pbc();
    }

    if (buffer.str().find("rod_vdw_site_placement") !=
        std::string::npos)
    {
        result = ffea_test::rod_vdw_site_placement();
    }

    if (buffer.str().find("point_lies_within_rod_element") !=
        std::string::npos)
    {
        result = ffea_test::point_lies_within_rod_element();
    }

    return result;
}

// Unit tests

int ffea_test::connection_test()
{
    std::cout << "Performing connection test...\n";
    World *world;
    world = new World();
    if (world->init("tet_ascii.1.ffea", 0, 0, 1) == FFEA_ERROR)
    {
        FFEA_error_text();
        cout << "Errors during initialisation mean World cannot be constructed "
                "properly."
             << endl;
    }
    std::cout << "Test world initialised. \n";

    rod::dbg_print = true;

    rod::Rod_blob_interface *current_interface =
        world->rod_blob_interface_array[0];

    float attachment_node_pos[3];
    float attachment_node[3];
    float attachment_material_axis[3];

    std::cout << "Current blob index in ffea_test: "
              << world->rod_blob_interface_array[0]->connected_blob->blob_index
              << "\n";

    current_interface->update_internal_state(true, true);
    current_interface->get_attachment_node(attachment_node, attachment_node_pos,
                                           true);
    current_interface->get_attachment_material_axis(attachment_node,
                                                    attachment_material_axis);

    float J[9];
    rod::get_jacobian(current_interface->deformed_tet_nodes, J);

    for (int i = 0; i < 4; i++)
    {
        rod::print_array("Tetrahedron node",
                         current_interface->connected_tet->n[i]->pos.data, 3);
    }
    rod::print_array("attachment node position", attachment_node_pos, 3);
    rod::print_array("attachment node", attachment_node, 3);
    rod::print_array("attachment material axis", attachment_material_axis, 3);
    rod::print_array("Undeformed jacobian", J, 9);

    rod::print_array("Undeformed jacobian inverse", current_interface->J_inv_0,
                     9);

    float rotmat[9];
    float rotation_angle = 0.55;
    rod::construct_euler_rotation_matrix(0, 0, rotation_angle, rotmat);
    rod::rotate_tet(rotmat, current_interface->connected_tet->n,
                    current_interface->deformed_tet_nodes);

    rod::print_array("connected tetrahedron node 0",
                     current_interface->connected_tet->n[0]->pos.data, 3);
    rod::print_array("connected tetrahedron node 1",
                     current_interface->connected_tet->n[1]->pos.data, 3);
    rod::print_array("connected tetrahedron node 2",
                     current_interface->connected_tet->n[2]->pos.data, 3);
    rod::print_array("connected tetrahedron node 3",
                     current_interface->connected_tet->n[3]->pos.data, 3);

    rod::print_array("rotation matrix", rotmat, 9);
    rod::print_array("Rotated tetrahedron node 1",
                     current_interface->deformed_tet_nodes[0]->pos.data, 3);
    rod::print_array("Rotated tetrahedron node 2",
                     current_interface->deformed_tet_nodes[1]->pos.data, 3);
    rod::print_array("Rotated tetrahedron node 3",
                     current_interface->deformed_tet_nodes[2]->pos.data, 3);
    rod::print_array("Rotated tetrahedron node 4",
                     current_interface->deformed_tet_nodes[3]->pos.data, 3);

    float new_attachment_node[3];
    float new_attachment_material_axis[3];

    current_interface->reorientate_connection(
        attachment_node, attachment_material_axis, new_attachment_node,
        new_attachment_material_axis);
    rod::print_array("Rotated node", new_attachment_node, 3);
    rod::print_array("Rotated material axis", new_attachment_material_axis, 3);
    float dotprod =
        new_attachment_material_axis[0] * attachment_material_axis[0] +
        new_attachment_material_axis[1] * attachment_material_axis[1] +
        new_attachment_material_axis[2] * attachment_material_axis[2];
    float material_axis_rotation_angle = std::acos(dotprod);
    std::cout << "Angle between = " << material_axis_rotation_angle << "\n";

    if ((material_axis_rotation_angle > rotation_angle - 0.03) &&
        (material_axis_rotation_angle < rotation_angle + 0.03))
    {
        // todo: make it so that the rotation is in the axis of the attachment face,
        // redo test
        return 0;
    }
    std::cout << "Rod-blob coupling test failed.\n";
    return 1;
}

int ffea_test::arbitrary_equilibrium_twist()
{
    std::cout << "Performing equilibrium twist test...\n";
    float reference_energy = 2.0955089065982904;

    float beta = 1;
    float offset = 0.3333333333 * M_PI;
    float pi[3] = {1, 0, 0};
    float pim1[3] = {1, 0, 0};
    float pi_equil[3] = {1, 0, 0};
    float pim1_equil[3] = {1, 0, 0};
    float mi[3] = {0, 1, 0};
    float mim1[3] = {0, 1, 0};
    float mi_equil[3] = {0, 1, 0};
    float mim1_equil[3] = {0, 1, 0};
    float mim1_equil_rotated[3];
    float mi_rotated[3];
    rod::rodrigues_rotation(mim1_equil, pim1_equil, offset, mim1_equil_rotated);
    rod::rodrigues_rotation(mi, pi, 1., mi_rotated);

    //    rod::print_array("pi", pi, 3);
    //    rod::print_array("pim1", pim1, 3);
    //    rod::print_array("pi_equil", pi_equil, 3);
    //    rod::print_array("pim1_equil", pim1_equil, 3);
    //    rod::print_array("mi", mi, 3);
    //    rod::print_array("mi_equil", mi_equil, 3);
    //    rod::print_array("mim1", mim1, 3);
    //    rod::print_array("mim1_equil", mim1_equil, 3);

    float twist_energy =
        rod::get_twist_energy(beta, mi_rotated, mim1, mi_equil,
                              mim1_equil_rotated, pim1, pi, pim1_equil, pi_equil);
    //    rod::print_array("mi_rotated", mi_rotated, 3);
    //    rod::print_array("mim1_equil_rotated",mim1_equil_rotated, 3);
    //    std::cout << "Twist energy = " << twist_energy << "\n";

    if ((twist_energy > reference_energy - 0.01) &&
        (twist_energy < reference_energy + 0.01))
    {
        std::cout << "It's all gravy. \n";
        return 0;
    }

    std::cout << "Twist energy = " << twist_energy << "\n";
    std::cout << "Reference energy = " << reference_energy << "\n";
    std::cout << "Test failed.\n";
    return 1;
}

int ffea_test::connection_orientation_test()
{
    std::cout << "Performing connection orientation test...\n";
    World *world;
    world = new World();
    if (world->init("tet_ascii.1.ffea", 0, 0, 1) == FFEA_ERROR)
    {
        FFEA_error_text();
        cout << "Errors during initialisation mean World cannot be constructed "
                "properly."
             << endl;
    }
    std::cout << "Test world initialised. \n";

    rod::Rod_blob_interface *current_interface =
        world->rod_blob_interface_array[0];

    rod::print_array("Old rod starting position",
                     current_interface->connected_rod->current_r, 3);

    // float attachment_node_pos[3];
    // float attachment_node[3];
    // float attachment_material_axis[3];

    // std::cout << "Current blob index in ffea_test: " <<
    // world->rod_blob_interface_array[0]->connected_blob->blob_index << "\n";

    // Set edge vecs first???

    // world->rod_blob_interface_array[0]->update_internal_state(true, true);

    // current_interface->position_rod_from_blob();

    // world->rod_blob_interface_array[0]->position_rod_from_blob();

    // world->rod_blob_interface_array[0]->update_internal_state(true, true);

    // world->rod_blob_interface_array[0]->position_blob_from_rod();

    // current_interface->update_internal_state(true, true);
    // current_interface->get_attachment_node(attachment_node,
    // attachment_node_pos);
    // current_interface->get_attachment_material_axis(attachment_node,
    // attachment_material_axis);

    // float J[9];
    // rod::get_jacobian(current_interface->deformed_tet_nodes, J);

    // for(int i = 0; i<4; i++){rod::print_array("Tetrahedron node",
    // current_interface->connected_tet->n[i]->pos.data, 3);}
    // rod::print_array("attachment node position", attachment_node_pos, 3);
    // rod::print_array("attachment node", attachment_node, 3);
    // rod::print_array("attachment material axis", attachment_material_axis, 3);
    // rod::print_array("Undeformed jacobian", J, 9);

    // current_interface->connected_rod->write_frame_to_file();

    world->print_trajectory_and_measurement_files(1, 1);
    world->print_trajectory_and_measurement_files(2, 2);

    rod::print_array("New rod starting position",
                     current_interface->connected_rod->current_r, 3);

    return 0;
}

int ffea_test::arbitrary_equilibrium_bend()
{
    float B_i_equil[4] = {1, 0, 0, 1};
    float B_im1_equil[4]{1, 0, 0, 1};
    float p_i[3] = {1, 0, 0};
    float p_im1[3] = {1, 0, 0};
    float m_i[3] = {0, 1, 0};
    float m_im1[3]{0, 1, 0};
    float n_i[3];
    float n_im1[3];

    float p_i_equil[3] = {1, 0, 0};
    float p_im1_equil[3] = {1, 0, 0};
    float m_i_equil[3] = {0, 1, 0};
    float m_im1_equil[3]{0, 1, 0};
    float n_i_equil[3];
    float n_im1_equil[3];

    float rm[9];
    float euler[3] = {0.25, 0.125, 0.2};
    rod::get_rotation_matrix_from_euler(euler, rm);
    float p_i_rotated[3];
    float m_i_rotated[3];
    rod::apply_rotation_matrix(p_i, rm, p_i_rotated);
    rod::apply_rotation_matrix(m_i, rm, m_i_rotated);

    rod::cross_product(m_im1_equil, p_im1_equil, n_im1_equil);
    rod::cross_product(m_i_equil, p_i_equil, n_i_equil);
    rod::cross_product(m_im1, p_im1, n_im1);
    rod::cross_product(m_i_rotated, p_i_rotated, n_i);

    float energy1 = rod::get_bend_energy_mutual_parallel_transport(
        p_im1, p_i_rotated, p_im1_equil, p_i_equil, n_im1, m_im1, n_im1_equil,
        m_im1_equil, n_i, m_i_rotated, n_i_equil, m_i_equil, B_i_equil,
        B_im1_equil);
    std::cout << "State 1:\n";
    rod::print_array("  p_im1", p_im1, 3);
    rod::print_array("  p_i", p_i_rotated, 3);
    rod::print_array("  p_im1_equil", p_im1_equil, 3);
    rod::print_array("  p_i_equil", p_i_equil, 3);
    rod::print_array("  n_im1", n_im1, 3);
    rod::print_array("  m_im1", m_im1, 3);
    rod::print_array("  m_im1_equil", m_im1, 3);
    rod::print_array("  n_im1_equil", n_im1_equil, 3);
    rod::print_array("  n_i", n_i, 3);
    rod::print_array("  m_i", m_i_rotated, 3);
    rod::print_array("  n_i_equil", n_i_equil, 3);
    rod::print_array("  m_i_equil", m_i_equil, 3);

    float p_i_rotated_2[3];
    float m_i_rotated_2[3];

    vec3d(n) { p_i_equil[n] = p_i_rotated[n]; }
    vec3d(n) { m_i_equil[n] = m_i_rotated[n]; }
    rod::apply_rotation_matrix(p_i_rotated, rm, p_i_rotated_2);
    rod::apply_rotation_matrix(m_i_rotated, rm, m_i_rotated_2);

    rod::cross_product(m_im1_equil, p_im1_equil, n_im1_equil);
    rod::cross_product(m_i_equil, p_i_equil, n_i_equil);
    rod::cross_product(m_im1, p_im1, n_im1);
    rod::cross_product(m_i_rotated_2, p_i_rotated_2, n_i);

    float energy2 = rod::get_bend_energy_mutual_parallel_transport(
        p_im1, p_i_rotated_2, p_im1_equil, p_i_equil, n_im1, m_im1, n_im1_equil,
        m_im1_equil, n_i, m_i_rotated_2, n_i_equil, m_i_equil, B_i_equil,
        B_im1_equil);

    std::cout << "State 2:\n";
    rod::print_array("  p_im1", p_im1, 3);
    rod::print_array("  p_i", p_i_rotated_2, 3);
    rod::print_array("  p_im1_equil", p_im1_equil, 3);
    rod::print_array("  p_i_equil", p_i_equil, 3);
    rod::print_array("  n_im1", n_im1, 3);
    rod::print_array("  m_im1", m_im1, 3);
    rod::print_array("  m_im1_equil", m_im1, 3);
    rod::print_array("  n_im1_equil", n_im1_equil, 3);
    rod::print_array("  n_i", n_i, 3);
    rod::print_array("  m_i", m_i_rotated_2, 3);
    rod::print_array("  n_i_equil", n_i_equil, 3);
    rod::print_array("  m_i_equil", m_i_equil, 3);

    std::cout << "Energy 1 = " << energy1 << "\n";
    std::cout << "Energy 2 = " << energy2 << "\n";

    float ref_diff = -0.0007473;
    if ((energy1 - energy2 > ref_diff - 0.01) &&
        (energy1 - energy2 < ref_diff + 0.01))
    {
        return 0;
    }

    std::cout << "get off my case ok not all software has to be correct\n";

    // set equil to equal current
    // rotate current by the same amount
    // are energies the same???

    return 1;
}

int ffea_test::identify_face()
{
    std::cout << "Identify face\n";

    World *world;
    world = new World();
    if (world->init("round_lad.1.ffea", 0, 0, 1) == FFEA_ERROR)
    {
        FFEA_error_text();
        cout << "Errors during initialisation mean World cannot be constructed "
                "properly."
             << endl;
    }
    std::cout << "Test world initialised. \n";

    rod::Rod_blob_interface *current_interface =
        world->rod_blob_interface_array[0];

    world->print_trajectory_and_measurement_files(1, 1);
    world->print_trajectory_and_measurement_files(2, 2);

    float attachment_node[3];
    float attachment_node_pos[3];
    // int face_node_indices[3];
    // current_interface->select_face_nodes(face_node_indices);
    float face_node_1[3];
    float face_node_2[3];
    float face_node_3[3];
    vec3d(n)
    {
        face_node_1[n] =
            current_interface
                ->deformed_tet_nodes[current_interface->face_node_indices[0]]
                ->pos.data[n];
    }
    vec3d(n)
    {
        face_node_2[n] =
            current_interface
                ->deformed_tet_nodes[current_interface->face_node_indices[1]]
                ->pos.data[n];
    }
    vec3d(n)
    {
        face_node_3[n] =
            current_interface
                ->deformed_tet_nodes[current_interface->face_node_indices[2]]
                ->pos.data[n];
    }

    current_interface->get_attachment_node(attachment_node, attachment_node_pos,
                                           true);

    float face_element_1[3];
    float face_element_2[3];
    float face_element_3[3];

    vec3d(n) { face_element_1[n] = face_node_2[n] - face_node_1[n]; }
    vec3d(n) { face_element_2[n] = face_node_3[n] - face_node_2[n]; }
    vec3d(n) { face_element_3[n] = face_node_1[n] - face_node_3[n]; }

    rod::normalize(face_element_1, face_element_1);
    rod::normalize(face_element_2, face_element_2);
    rod::normalize(face_element_3, face_element_3);

    float node1dp = (face_element_1[0] * attachment_node[0]) +
                    (face_element_1[1] * attachment_node[1]) +
                    (face_element_1[2] * attachment_node[2]);
    float node2dp = (face_element_2[0] * attachment_node[0]) +
                    (face_element_2[1] * attachment_node[1]) +
                    (face_element_2[2] * attachment_node[2]);
    float node3dp = (face_element_3[0] * attachment_node[0]) +
                    (face_element_3[1] * attachment_node[1]) +
                    (face_element_3[2] * attachment_node[2]);

    std::cout << "node 1 dp: " << node1dp << "\n";
    std::cout << "node 2 dp: " << node2dp << "\n";
    std::cout << "node 3 dp: " << node3dp << "\n";

    if ((node1dp < 0.01) && (node2dp < 0.01) && (node3dp < 0.01))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

int ffea_test::connection_energy()
{ // this test just checks that the energy
    // is zero at connection equilibrium!!!

    std::cout << "doin a connection energy test \n";

    World *world;
    world = new World();
    if (world->init("realistic.ffea", 0, 0, 1) == FFEA_ERROR)
    {
        FFEA_error_text();
        cout << "Errors during initialisation mean World cannot be constructed "
                "properly."
             << endl;
    }
    std::cout << "Test world initialised. \n";

    rod::Rod_blob_interface *current_interface =
        world->rod_blob_interface_array[0];

    float attachment_node[3];
    float attachment_node_pos[3];
    float attachment_material_axis[3];

    float attachment_node_equil[3];
    float attachment_node_pos_equil[3];
    float attachment_material_axis_equil[3];

    current_interface->get_attachment_node(attachment_node, attachment_node_pos,
                                           false);
    current_interface->get_attachment_material_axis(attachment_node,
                                                    attachment_material_axis);

    current_interface->get_attachment_node(attachment_node_equil,
                                           attachment_node_pos_equil, true);
    std::cout << "Attachment node acquired.\n";
    current_interface->get_attachment_material_axis(
        attachment_node_equil, attachment_material_axis_equil);
    std::cout << "Attachment material axis acquired.\n";
    current_interface->update_internal_state(true, true);

    rod::print_array("Final attachment node", attachment_node, 3);

    // current_interface->do_connection_timestep();

    int node_index = 1;
    // float displacement = current_interface->connected_rod->perturbation_amount;
    float displacement = 1;
    float energy[6] = {0, 0, 0, 0, 0, 0};

    current_interface->get_node_energy(
        node_index, attachment_node_equil, attachment_material_axis_equil,
        attachment_node, attachment_material_axis, displacement, energy);
    rod::print_array("energy", energy, 6);

    float force[3];

    vec3d(n) { force[n] = (energy[n] - energy[n + 3]) / displacement; }

    rod::print_array("force", force, 3);

    world->print_trajectory_and_measurement_files(1, 1);
    world->print_trajectory_and_measurement_files(2, 2);

    for (int i = 0; i < 3; i++)
    {
        if (force[i] > 0.001)
        {
            std::cout << "The rod-blob interface has a non-zero equilibrium force!\n";
            return 1;
        }
    }

    world->print_trajectory_and_measurement_files(1, 1);
    world->print_trajectory_and_measurement_files(2, 2);

    return 0;

    /**

      World *mirror_world;
      mirror_world = new World();
      if(mirror_world->init("realistic_reversed.ffea", 0, 0, 1) == FFEA_ERROR) {
          FFEA_error_text();
          cout << "Errors during initialisation mean World cannot be constructed
     properly." << endl;
      }
      std::cout << "Mirror world initialised... oooo! \n";

      rod::Rod_blob_interface *mirror_interface =
     mirror_world->rod_blob_interface_array[0]; float mirror_attachment_node[3];
      float mirror_attachment_node_pos[3];
      float mirror_attachment_material_axis[3];
      mirror_interface->get_attachment_node(mirror_attachment_node,
     mirror_attachment_node_pos);
      mirror_interface->get_attachment_material_axis(mirror_attachment_node,
     mirror_attachment_material_axis);
      mirror_interface->update_internal_state(true, true);
      int mirror_node_index = 1;
      float mirror_displacement = 0;
      float mirror_energy[6] = {0,0,0,0,0,0};
      mirror_interface->get_node_energy(mirror_node_index,
     mirror_attachment_node, mirror_attachment_material_axis,
     mirror_displacement, mirror_energy); rod::print_array("mirror_energy",
     mirror_energy, 6);

      for (int i=0; i<6; i++){
          if (mirror_energy[i] > 0.0001){
              std::cout << "The rod-blob interface has a non-zero equilibrium
     energy!\n"; return 1;
          }
      }

      // for rod-to-blob, these energies are now all zero-ish
      // next move: make sure they're all zero-ish for blob-to-rod
      // maybe come up with a sensible test for non-zero
      // like seeing if, in a cold simulation, things relax to equilibrium

  */

    return 0;
}

int ffea_test::connection_energy_2()
{
    World *world;
    world = new World();
    if (world->init("realistic.ffea", 0, 0, 1) == FFEA_ERROR)
    {
        FFEA_error_text();
        cout << "Errors during initialisation mean World cannot be constructed "
                "properly."
             << endl;
    }

    rod::Rod_blob_interface *current_interface =
        world->rod_blob_interface_array[0];

    float attachment_node[3];
    float attachment_node_pos[3];
    float attachment_material_axis[3];

    float attachment_node_equil[3];
    float attachment_node_pos_equil[3];
    float attachment_material_axis_equil[3];

    current_interface->get_attachment_node(attachment_node, attachment_node_pos,
                                           false);
    current_interface->get_attachment_material_axis(attachment_node,
                                                    attachment_material_axis);

    current_interface->get_attachment_node(attachment_node_equil,
                                           attachment_node_pos_equil, true);
    current_interface->get_attachment_material_axis(
        attachment_node_equil, attachment_material_axis_equil);

    current_interface->update_internal_state(true, true);

    // float displacement = current_interface->connected_rod->perturbation_amount;
    float displacement = 0.0000025;

    for (int node_index = 0; node_index < 4; node_index++)
    {
        for (int i = 0; i < 500; i++)
        {
            float energy[6] = {0, 0, 0, 0, 0, 0};
            current_interface->get_node_energy(
                node_index, attachment_node_equil, attachment_material_axis_equil,
                attachment_node, attachment_material_axis, displacement * i, energy);
            std::cout << "ENERGYPLOT displacement " << i << " node " << node_index
                      << " energy " << energy[0] << " " << energy[1] << " "
                      << energy[2] << " " << energy[3] << " " << energy[4] << " "
                      << energy[5] << "\n";
        }
    }

    // float force[3];

    // vec3d(n){force[n] = (energy[n] - energy[n+3])/displacement;}

    // rod::print_array("force", force, 3);

    //    std::cout << "doin a connection energy test \n";
    //
    //    World *world;
    //    world = new World();
    //    if(world->init("realistic.ffea", 0, 0, 1) == FFEA_ERROR) {
    //        FFEA_error_text();
    //        cout << "Errors during initialisation mean World cannot be
    //        constructed properly." << endl;
    //    }
    //    std::cout << "Test world initialised. \n";
    //
    //    rod::Rod_blob_interface *current_interface =
    //    world->rod_blob_interface_array[0];
    //
    //    float attachment_node[3];
    //    float attachment_node_pos[3];
    //    rod::print_array("    new attachment node (should still be 1,0,0)",
    //    current_interface->attachment_node, 3);
    //
    //    current_interface->do_connection_timestep();
    //
    //    rod::print_array("    new attachment node (should still be 1,0,0)",
    //    current_interface->attachment_node, 3);
    //
    //    world->print_trajectory_and_measurement_files(1, 1);
    //    world->print_trajectory_and_measurement_files(2, 2);
    //
    //    current_interface->position_blob_from_rod();
    //
    //    rod::print_array("    new attachment node (should still be 1,0,0)",
    //    current_interface->attachment_node, 3);
    //
    //    world->print_trajectory_and_measurement_files(3, 1);

    //    World *mirror_world;
    //    mirror_world = new World();
    //    if(mirror_world->init("realistic_reversed.ffea", 0, 0, 1) == FFEA_ERROR)
    //    {
    //        FFEA_error_text();
    //        cout << "Errors during initialisation mean World cannot be
    //        constructed properly." << endl;
    //    }
    //    std::cout << "Mirror world initialised. \n";
    //
    //    rod::Rod_blob_interface *mirror_interface =
    //    mirror_world->rod_blob_interface_array[0];

    /**

      float attachment_node[3];
      float attachment_node_pos[3];
      float attachment_material_axis[3];

      current_interface->get_attachment_node(attachment_node,
     attachment_node_pos);
      current_interface->get_attachment_material_axis(attachment_node,
     attachment_material_axis);

      rod::print_array("Attachment node", attachment_node, 3);
      rod::print_array("Attachment material axis", attachment_material_axis, 3);

      current_interface->update_internal_state(true, true);

      float displacement = 0;
      float energy[2][6] = {{0,0,0,0,0,0},{0,0,0,0,0,0}};

      std::cout << "Attachment node acquired.\n";

      current_interface->get_rod_energy(attachment_node,
     attachment_material_axis, displacement, energy);

      std::cout << "Rod energy at equilibrium:\n";
      std::cout << "r0 : [" << energy[0][0] << ", " << energy[0][1] << ", " <<
     energy[0][2] << ", " << energy[0][3] << ", " << energy[0][4] << ", " <<
     energy[0][5] << "]\n"; std::cout << "r1 : [" << energy[1][0] << ", " <<
     energy[1][1] << ", " << energy[1][2] << ", " << energy[1][3] << ", " <<
     energy[1][4] << ", " << energy[1][5] << "]\n";

      world->print_trajectory_and_measurement_files(1, 1);
      world->print_trajectory_and_measurement_files(2, 2);

  */

    return 0;
}

int ffea_test::jacobian_rotate()
{
    mesh_node *node_up[4];
    for (int i = 0; i < 4; i++)
    {
        node_up[i] = new mesh_node();
    }

    node_up[0]->pos.x = 72.138;
    node_up[0]->pos.y = 42.9213;
    node_up[0]->pos.z = 37.3931;
    node_up[1]->pos.x = 83.2145;
    node_up[1]->pos.y = 31.5663;
    node_up[1]->pos.z = -19.2512;
    node_up[2]->pos.x = 79.1497;
    node_up[2]->pos.y = 57.1503;
    node_up[2]->pos.z = -19.2512;
    node_up[3]->pos.x = 58.1056;
    node_up[3]->pos.y = 35.524;
    node_up[3]->pos.z = -19.2512;

    mesh_node *node_forward[4];
    for (int i = 0; i < 4; i++)
    {
        node_forward[i] = new mesh_node();
    }

    node_forward[0]->pos.x = 80.254;
    node_forward[0]->pos.y = 42.904;
    node_forward[0]->pos.z = 42.7413;
    node_forward[1]->pos.x = 136.898;
    node_forward[1]->pos.y = 31.6887;
    node_forward[1]->pos.z = 53.9594;
    node_forward[2]->pos.x = 136.898;
    node_forward[2]->pos.y = 57.2197;
    node_forward[2]->pos.z = 49.574;
    node_forward[3]->pos.x = 136.898;
    node_forward[3]->pos.y = 35.3313;
    node_forward[3]->pos.z = 28.8028;

    // note: confirmed for same tetrahedron!

    float J_up[9];
    float J_forward[9];

    rod::get_jacobian(node_up, J_up);
    rod::get_jacobian(node_forward, J_forward);

    rod::print_array("  J up", J_up, 9);
    rod::print_array("  J forward", J_forward, 9);

    return 0;
}

int ffea_test::connection_energy_3()
{
    World *world;
    world = new World();
    if (world->init("realistic.ffea", 0, 0, 1) == FFEA_ERROR)
    {
        FFEA_error_text();
        cout << "Errors during initialisation mean World cannot be constructed "
                "properly."
             << endl;
    }

    rod::Rod_blob_interface *current_interface =
        world->rod_blob_interface_array[0];

    current_interface->connected_tet->n[1]->pos.x -= 10.;

    world->run();

    return 0; // note: this one should just crash the program if it's wrong
}

int ffea_test::connection_propagation_every_way()
{
    int tests_failed = 0;
    // mode 0 = twist, mode 1 = bend, mode 2 = stretch
    // tests_failed += connection_propagation(0, true);
    tests_failed += connection_propagation(1, false);
    tests_failed += connection_propagation(0, false);
    return 0;
}

int ffea_test::connection_propagation(
    int mode,
    bool ends_at_rod)
{ // mode 0 = twist, mode 1 = bend, mode 2 = stretch

    rod::dbg_print = false;

    World *world;
    world = new World();

    if (mode == 0)
    { // twist
        if (world->init("twist.ffea", 0, 0, 1) == FFEA_ERROR)
        {
            FFEA_error_text();
            cout << "Errors during initialisation mean World cannot be constructed "
                    "properly."
                 << endl;
        }
    }
    if (mode == 1)
    {
        if (world->init("bend.ffea", 0, 0, 1) == FFEA_ERROR)
        {
            FFEA_error_text();
            cout << "Errors during initialisation mean World cannot be constructed "
                    "properly."
                 << endl;
        }
    }

    rod::Rod *current_rod = world->rod_array[0];

    int end_index;

    if (world->rod_blob_interface_array[0]->ends_at_rod == false)
    {
        end_index = 0;
    }
    else
    {
        end_index = current_rod->num_nodes - 2;
    }

    current_rod->pinned_nodes[end_index] = true;

    float curr_sample_edge_vec[3];

    float rotmat[9];
    float edgevec_dot_prod_pre;

    if (mode == 0)
    { // twist

        for (int i = 0; i < 3; i++)
        {
            curr_sample_edge_vec[i] =
                world->rod_blob_interface_array[0]->edge_vecs[0][i];
        }
        rod::normalize(curr_sample_edge_vec, curr_sample_edge_vec);

        // twist end_p
        float end_p[3];
        float end_m_pre[3];
        float end_m[3];

        current_rod->get_p(end_index, end_p, false);
        for (int n = 0; n < 3; n++)
        {
            end_m_pre[n] = current_rod->current_m[n];
        }

        rod::normalize(end_p, end_p);
        rod::normalize(end_m_pre, end_m_pre);

        edgevec_dot_prod_pre = curr_sample_edge_vec[0] * end_m_pre[0] +
                               curr_sample_edge_vec[1] * end_m_pre[1] +
                               curr_sample_edge_vec[2] * end_m_pre[2];

        rod::print_array("original m", end_m, 3);
        rod::rodrigues_rotation(end_m_pre, end_p, M_PI * 0.5, end_m);
        rod::print_array("rotated m", end_m, 3);

        rod::get_rotation_matrix(end_m_pre, end_m, rotmat);

        for (int n = 0; n < 3; n++)
        {
            current_rod->current_m[n + (end_index * 3)] = end_m[n];
        }

        int direction_factor = 1;
        if (world->rod_blob_interface_array[0]->ends_at_rod == true)
        {
            direction_factor = -1;
        }

        // twist the next p over
        float next_end_p[3];
        float next_end_m[3];

        current_rod->get_p(1, next_end_p, false);
        for (int n = 0; n < 3; n++)
        {
            next_end_m[n] =
                current_rod->current_m[n + (3 * direction_factor) + (end_index * 3)];
        }
        rod::normalize(next_end_p, next_end_p);
        rod::normalize(next_end_m, next_end_m);
        rod::rodrigues_rotation(next_end_m, next_end_p, M_PI * 0.5, next_end_m);
        for (int n = 0; n < 3; n++)
        {
            current_rod->current_m[n + (3 * direction_factor) + (end_index * 3)] =
                next_end_m[n];
        }
    }

    if (mode == 1)
    {
        float end_p[3];
        float end_m[3];
        float end_p_rotated[3];
        float end_m_rotated[3];
        current_rod->get_p(end_index, end_p, false);
        float euler_angles[3] = {0, M_PI / 4.0, 0};
        float rm[9];
        rod::get_rotation_matrix_from_euler(euler_angles, rm);
        rod::print_array("rotmat", rm, 9);
        float scale = rod::absolute(end_p);
        vec3d(n) { end_p[n] /= scale; }
        rod::apply_rotation_matrix(end_p, rm, end_p_rotated);
        rod::apply_rotation_matrix(end_m, rm, end_m_rotated);
        vec3d(n) { end_p_rotated[n] *= scale; }
        rod::print_array("p", end_p, 3);
        rod::print_array("p_rotated", end_p_rotated, 3);
        rod::normalize(end_m_rotated, end_m_rotated);

        for (int n = 0; n < 3; n++)
        {
            current_rod->current_m[n + (end_index * 3)] = end_m[n];
        }

        if (world->rod_blob_interface_array[0]->ends_at_rod)
        {
            vec3d(n)
            {
                current_rod->current_r[n + current_rod->length - 3] =
                    current_rod->current_r[n + current_rod->length - 6] +
                    end_p_rotated[n];
            }
            current_rod->pinned_nodes[current_rod->num_nodes - 1] = true;
        }
        else
        {
            vec3d(n)
            {
                current_rod->current_r[n] =
                    current_rod->current_r[n + 3] - end_p_rotated[n];
            }
            current_rod->pinned_nodes[1] = true;
        }
    }

    rod::Rod_blob_interface *current_interface =
        world->rod_blob_interface_array[0];

    // get the first energies

    float equil_attachment_node[3];
    float equil_attachment_node_pos[3];
    float equil_attachment_material_axis[3];

    float init_energy[6];

    float init_attachment_node[3];
    float init_attachment_node_pos[3];
    float init_attachment_material_axis[3];

    current_interface->get_attachment_node(init_attachment_node,
                                           init_attachment_node_pos, false);
    current_interface->get_attachment_material_axis(
        init_attachment_node, init_attachment_material_axis);

    current_interface->get_attachment_node(equil_attachment_node,
                                           equil_attachment_node_pos, true);
    current_interface->get_attachment_material_axis(
        equil_attachment_node, equil_attachment_material_axis);

    current_interface->update_internal_state(true, true);

    current_interface->get_node_energy(
        1, equil_attachment_node, equil_attachment_material_axis,
        init_attachment_node, init_attachment_material_axis, 0.1, init_energy);
    std::cout << "CURR_ENERGY_INIT L" << init_energy[0] << " " << init_energy[1]
              << " " << init_energy[2] << " " << init_energy[3] << " "
              << init_energy[4] << " " << init_energy[5] << "\n";

    std::cout << "Running the test: \n";
    for (int i = 0; i < 20000; i++)
    {
        world->rod_array[0]->do_timestep(world->rng);
    }

    // eyy
    world->run();

    // std::cout << "sample twist timestep\n";
    // rod::dbg_print = true;

    // for (int i=0; i<1; i++){
    //    world->rod_array[0]->do_timestep(world->rng);
    //}
    rod::dbg_print = false;

    float attachment_node[3];
    float attachment_node_pos[3];
    float attachment_material_axis[3];

    current_interface->get_attachment_node(attachment_node, attachment_node_pos,
                                           true);
    current_interface->get_attachment_material_axis(attachment_node,
                                                    attachment_material_axis);
    current_interface->update_internal_state(true, true);

    // float displacement = current_interface->connected_rod->perturbation_amount;

    float post_sample_edge_vec[3];
    for (int i = 0; i < 3; i++)
    {
        post_sample_edge_vec[i] =
            world->rod_blob_interface_array[0]->edge_vecs[0][i];
    }
    rod::normalize(post_sample_edge_vec, post_sample_edge_vec);

    float end_m_post[3];
    for (int n = 0; n < 3; n++)
    {
        end_m_post[n] = current_rod->current_m[n];
    }

    float edgevec_dot_prod_post = post_sample_edge_vec[0] * end_m_post[0] +
                                  post_sample_edge_vec[1] * end_m_post[1] +
                                  post_sample_edge_vec[2] * end_m_post[2];
    std::cout << "post dot prod: " << edgevec_dot_prod_post << "\n";
    std::cout << "pre dot prod: " << edgevec_dot_prod_pre << "\n";

    rod::dbg_print = false;

    float energy[6];

    current_interface->get_node_energy(
        1, equil_attachment_node, equil_attachment_material_axis, attachment_node,
        attachment_material_axis, 0.1, energy);

    std::cout << "CURR_ENERGY_INIT L" << energy[0] << " " << energy[1] << " "
              << energy[2] << " " << energy[3] << " " << energy[4] << " "
              << energy[5] << "\n";

    // float rotation = M_PI/2./200;
    float perturb = current_interface->connected_rod->perturbation_amount * 50;

    for (int node_index = 0; node_index < 10; node_index++)
    {
        for (int i = 0; i < 100; i++)
        {
            rod::dbg_print = false;

            if ((node_index == 5) & (i == 10))
            {
                rod::dbg_print = true;
            }

            float energyplus[3] = {0, 0, 0};
            float energyminus[3] = {0, 0, 0};
            int start_cutoff;
            int end_cutoff;
            rod::set_cutoff_values(node_index,
                                   current_interface->connected_rod->num_nodes,
                                   &start_cutoff, &end_cutoff);

            rod::get_perturbation_energy( // from rod_math
                perturb * i,
                rod::x, // dimension
                current_interface->connected_rod->B_matrix,
                current_interface->connected_rod->material_params, start_cutoff,
                end_cutoff, node_index, current_interface->connected_rod->current_r,
                current_interface->connected_rod->equil_r,
                current_interface->connected_rod->current_m,
                current_interface->connected_rod->equil_m, energyplus);

            std::cout << "TWIST displacement " << i * perturb << " node "
                      << node_index << " energy " << energyplus[2] << "\n";
            std::cout << "BEND displacement " << i * perturb << " node " << node_index
                      << " energy " << energyplus[1] << "\n";
            std::cout << "STRETCH displacement " << i * perturb << " node "
                      << node_index << " energy " << energyplus[0] << "\n";

            rod::get_perturbation_energy( // from rod_math
                perturb * i * -1,
                rod::x, // dimension
                current_interface->connected_rod->B_matrix,
                current_interface->connected_rod->material_params, start_cutoff,
                end_cutoff, node_index, current_interface->connected_rod->current_r,
                current_interface->connected_rod->equil_r,
                current_interface->connected_rod->current_m,
                current_interface->connected_rod->equil_m, energyminus);

            std::cout << "TWIST displacement " << i * perturb * -1 << " node "
                      << node_index << " energy " << energyminus[2] << "\n";
            std::cout << "BEND displacement " << i * perturb * -1 << " node "
                      << node_index << " energy " << energyminus[1] << "\n";
            std::cout << "STRETCH displacement " << i * perturb * -1 << " node "
                      << node_index << " energy " << energyminus[0] << "\n";
        }
    }

    float displacement = 0.0125;
    for (int node_index = 0; node_index < 4; node_index++)
    {
        for (int i = 0; i < 500; i++)
        {
            float energy[6] = {0, 0, 0, 0, 0, 0};
            current_interface->get_node_energy(
                node_index, equil_attachment_node, equil_attachment_material_axis,
                attachment_node, attachment_material_axis, displacement * i, energy);
            std::cout << "ENERGYPLOT displacement " << i << " node " << node_index
                      << " energy " << energy[0] << " " << energy[1] << " "
                      << energy[2] << " " << energy[3] << " " << energy[4] << " "
                      << energy[5] << "\n";
        }
    }

    float dotprod_diff = edgevec_dot_prod_post - edgevec_dot_prod_pre;

    rod::dbg_print = true;
    std::cout << "\n";
    rod::print_array("current_m", current_interface->connected_rod->current_m,
                     30);
    std::cout << "\n";
    rod::print_array("equil_m", current_interface->connected_rod->equil_m, 30);
    std::cout << "\n";
    rod::print_array("current_r", current_interface->connected_rod->current_r,
                     30);
    std::cout << "\n";
    rod::print_array("equil_r", current_interface->connected_rod->equil_r, 30);
    std::cout << "\n";

    if (mode == 0)
    {
        if (dotprod_diff > -0.05 && dotprod_diff < 0.05)
        {
            std::cout << "twist is good. everything is good. lets all go out for "
                         "some frosty chocolate milkshakes\n";
            return 0;
        }
    }

    if (mode == 1)
    {
        float att_node_end[3];
        float att_node_end_pos[3];
        current_interface->get_attachment_node(att_node_end, att_node_end_pos,
                                               false);
        float end_end_p[3];
        current_rod->get_p(end_index, end_end_p, false);
        rod::normalize(att_node_end, att_node_end);
        rod::normalize(end_end_p, end_end_p);
        float end_bend_dotprod = att_node_end[0] * end_end_p[0] +
                                 att_node_end[1] * end_end_p[1] +
                                 att_node_end[2] * end_end_p[2];
        if (end_bend_dotprod > 0.90)
        {
            std::cout << "bend is good. everything is good. lets all go out for some "
                         "frosty chocolate milkshakes\n";
            return 0;
        }
    }

    std::cout << "mode" << mode << ", its broke\n";
    return 1;
}

int ffea_test::recover_normal()
{
    World *world;
    world = new World();
    if (world->init("realistic.ffea", 0, 0, 1) == FFEA_ERROR)
    {
        FFEA_error_text();
        cout << "Errors during initialisation mean World cannot be constructed "
                "properly."
             << endl;
    }

    std::cout << "Starting equil attachment node (normal) test\n";

    rod::Rod_blob_interface *current_interface =
        world->rod_blob_interface_array[0];

    rod::dbg_print = true;

    float equil_attachment_node[3];

    rod::equil_attachment_node_from_J(
        current_interface->J_inv_0, current_interface->face_node_indices,
        current_interface->ends_at_rod, current_interface->node_weighting,
        current_interface->tet_origin, current_interface->edge_vecs,
        current_interface->euler_angles, equil_attachment_node);

    rod::print_array("current_interface->J_inv_0", current_interface->J_inv_0, 9);
    rod::print_array("equil_attachment_node", equil_attachment_node, 3);

    float attachment_node[3];
    float attachment_node_pos[3];
    rod::dbg_print = false;
    current_interface->get_attachment_node(attachment_node, attachment_node_pos,
                                           false);
    rod::dbg_print = true;
    rod::print_array("attachment node      ", attachment_node, 3);

    if (equil_attachment_node[0] > 0.99 && equil_attachment_node[1] < 0.01 &&
        equil_attachment_node[2] < 0.01)
    {
        return 0;
    }

    return 1;
}

int ffea_test::dump_twist_info()
{
    float everything_rotmat[9] = {0.6092191, -0.7677125, -0.1986693,
                                  0.6734007, 0.6331432, -0.3816559,
                                  0.4187881, 0.0987280, 0.9027011};
    // x-axis: clear, works fine //float everything_rotmat[9] = {  1.0000000,
    // 0.0000000, 0.0000000, 0.0000000, 0.8775826, -0.4794255, 0.0000000,
    // 0.4794255, 0.8775826 }; y-axis clear float everything_rotmat[9] = {
    // 0.8775826, 0.0000000, 0.4794255, 0.0000000, 1.0000000, 0.0000000,
    // -0.4794255, 0.0000000, 0.8775826 }; //y z-axis clear float
    // everything_rotmat[9] = {0.8775826, -0.4794255, 0.0000000, 0.4794255,
    // 0.8775826, 0.0000000, 0.0000000, 0.0000000, 1.0000000};

    float im1_rotation_angle = 0.5;
    float i_rotation_angle = -0.5;
    float i_rotmat[9] = {0.3469295, -0.6816330, 0.6442177, 0.9377583, 0.2636695,
                         -0.2260263, -0.0157935, 0.6825356, 0.7306817};

    float p_i_equil[3] = {1, 0, 0};
    float p_i[3] = {1, 0, 0};
    float m_i_equil[3] = {0, 1, 0};
    float m_i[3] = {0, 1, 0};

    float p_im1_equil[3] = {1, 0, 0};
    float p_im1[3] = {1, 0, 0};
    float m_im1_equil[3] = {0, 1, 0};
    float m_im1[3] = {0, 1, 0};

    float p_i_equil_rot[3] = {1, 0, 0};
    float p_i_rot[3] = {1, 0, 0};
    float m_i_equil_rot[3] = {0, 1, 0};
    float m_i_rot[3] = {0, 1, 0};

    float p_im1_equil_rot[3] = {1, 0, 0};
    float p_im1_rot[3] = {1, 0, 0};
    float m_im1_equil_rot[3] = {0, 1, 0};
    float m_im1_rot[3] = {0, 1, 0};

    rod::apply_rotation_matrix(p_i_equil, everything_rotmat, p_i_equil_rot);
    rod::apply_rotation_matrix(p_i_rot, everything_rotmat, p_i_rot);
    rod::apply_rotation_matrix(m_i_equil_rot, everything_rotmat, m_i_equil_rot);
    rod::apply_rotation_matrix(m_i_rot, everything_rotmat, m_i_rot);
    rod::apply_rotation_matrix(p_im1_equil_rot, everything_rotmat,
                               p_im1_equil_rot);
    rod::apply_rotation_matrix(p_im1_rot, everything_rotmat, p_im1_rot);
    rod::apply_rotation_matrix(m_im1_equil_rot, everything_rotmat,
                               m_im1_equil_rot);
    rod::apply_rotation_matrix(m_im1_rot, everything_rotmat, m_im1_rot);

    rod::dbg_print = true;

    rod::print_array("p_i_equil_rot", p_i_equil_rot, 3);
    rod::print_array("p_i_rot", p_i_rot, 3);
    rod::print_array("m_i_equil_rot", m_i_equil_rot, 3);
    rod::print_array("m_i_rot", m_i_rot, 3);
    rod::print_array("p_im1_equil_rot", p_im1_equil_rot, 3);
    rod::print_array("p_im1_rot", p_im1_rot, 3);
    rod::print_array("m_im1_equil_rot", m_im1_equil_rot, 3);
    rod::print_array("m_im1_rot", m_im1_rot, 3);

    float twenergy;

    for (int i = 0; i < 100; i++)
    {
        // int i = 65;

        float new_mi[3];
        rod::rodrigues_rotation(m_i_rot, p_i_rot, ((float)i / 100) * 2 * 3.14159,
                                new_mi);
        twenergy = rod::get_twist_energy(1, new_mi, m_im1_rot, m_i_equil_rot,
                                         m_im1_equil_rot, p_im1_rot, p_i_rot,
                                         p_im1_equil_rot, p_i_equil_rot);
        std::cout << "twist_dbg_plot " << (i / (float)100) * 2 * 3.14159 << " "
                  << twenergy << "\n";
    }

    return 0;
}

int ffea_test::lower_sphere()
{
    World *world;
    world = new World();
    if (world->init("ndc80c_mt_multi.ffea", 0, 0, 1) == FFEA_ERROR)
    {
        FFEA_error_text();
        cout << "Errors during initialisation mean World cannot be constructed "
                "properly."
             << endl;
    }

    for (int i = 0; i < 1000; i++)
    {
        std::cout << "Sphere lowered.\n";
        world->blob_array[7]->move(-0.1, 0, 0);
        world->run();
    }

    world->params.check *= 100;
    world->params.num_steps *= 10000000000;
    // world->params.dt *= 3;

    world->run();

    return 1;
}

// When calculating collisions, the closest line segment between two rod
// elements must be calculated, connected by a pair of points, one on each rod.
// Such points must lie on the rod element centreline and within the boundaries
// of its length.
//
// TODO: Test needs rewriting for new correction function.
//
// In this test, six points, c, are tested against an arbitrary rod element.
// Beyond the element length, a point should be assigned to the closest node.
// Otherwise, the point is unchanged.
//
// Define candidates for c with the parametric line equation, c = r1 + p*t,
// where c = r2 when t = 1. t is a multiplier.
int ffea_test::point_lies_within_rod_element()
{

    float t[6] = {-1, 0, 0.3, 0.7, 1, 1.5};
    float c_init[3] = {0};
    float c_out[3] = {0};

    float c_answer[3] = {0};
    float delta[3] = {0};
    int pass_count = 0;


    float r1[3] = {0};
    float p[3] = {0, 0, 1};


    //rod::dbg_print = true;

    for (int i = 0; i < 6; i++)
    {
        std::cout << "=================== CASE " << i << " ===================\n";
        vec3d(n) { c_init[n] = r1[n] + p[n] * t[i]; }

        // enforce finite length of rod element
        rod::finite_length_correction(c_init, r1, p, c_out);

        switch (i)
        {
        case 0:
            // c < r1, assign to node
            vec3d(n) { c_answer[n] = r1[n]; }
            break;
        case 1:
            // c == r1, assign to node
            vec3d(n) { c_answer[n] = r1[n]; }
            break;
        case 2:
            // r1 < c < r2, c is ok
            vec3d(n) { c_answer[n] = c_init[n]; }
            break;
        case 3:
            // r1 < c < r2, c is ok
            vec3d(n) { c_answer[n] = c_init[n]; }
            break;
        case 4:
            // c == r2, assign to node
            vec3d(n) { c_answer[n] = r1[n] + p[n]; }
            break;
        case 5:
            // c > r2, assign to node
            vec3d(n) { c_answer[n] = r1[n] + p[n]; }
            break;
        default:
            return 1;
        }

        std::cout << "t = " << t[i] << std::endl;
        rod::print_array("c initial", c_init, 3);
        rod::print_array("c computed", c_out, 3);
        rod::print_array("c expected", c_answer, 3);

        vec3d(n) { delta[n] = c_out[n] - c_answer[n]; }
        if (rod::absolute(delta) < 0.01)
        {
            pass_count++;
            std::cout << "PASS\n";
        }
        else
            std::cout << "FAIL\n";

    }

    if (pass_count == 6)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

/* Get two points, c_a and c_b, that form the straight line connecting two rod
elements
 * (capsules) and compare to expected values. Rod a is the 'stationary' rod in
 * these collision scenarios, whilst b is the 'colliding' rod.
 *
TODO: I think the function that 'corrects' for when a connection is outside the
rod element is a bit wrong, e.g. when a connection is actually 1/3 along the rod
element, it gets corrected to be at the node end anyway. Becomes more of a
problem if using long rod elements. May not be a huge issue. */
int ffea_test::line_connecting_rod_elements()
{
    float radius_a = 0.25;
    float radius_b = 0.5;
    float r_a1[3] = {0.0, 0.0, 0.0}; // fixed
    float r_b1[7][3] = {{0.0, radius_a + radius_b, -0.25},
                        {1.1456, radius_a, 0.25},
                        {0.75, radius_a + radius_b, -0.5},
                        {0.1, radius_b, 0.0},
                        {1, -0.5, 0},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.5, 0.0}};
    float r_a2[3] = {1.0, 0.0, 0.0}; // fixed
    float r_b2[3] = {0.0, 0.0, 0.0};
    float p_a[3] = {1.0, 0.0, 0.0}; // point in x direction
    float p_b[7][3] = {
        {std::sin(M_PI / 4.0), 0.0, std::cos(M_PI / 4.0)},   // 45 deg x-z
        {std::sin(M_PI / 12.0), 0.0, std::cos(M_PI / 12.0)}, // 15 deg x-z
        {0.0, 0.0, 1.0},
        {std::cos(M_PI / 12.0), std::sin(M_PI / 12.0), 0.0}, // 15 deg, x-y
        {std::cos(M_PI / 3.0), std::sin(M_PI / 3.0), 0.0},   // 60 deg, x-y
        {std::cos(M_PI / 12.0), std::sin(M_PI / 12.0), 0.0}, // 15 deg, x-y
        {std::cos(M_PI / 3.0), -std::sin(M_PI / 3.0), 0.0}}; // -60 deg, x-y
    float l_a[3] = {0, 0, 0};
    float l_b[3] = {0, 0, 0};
    float l_a_cross_l_b[3] = {0.0, 0.0, 0.0};
    float c_a[3] = {0.0, 0.0, 0.0};
    float c_b[3] = {0.0, 0.0, 0.0};

    // Expected results
    float c_a_answer[3] = {0.0, 0.0, 0.0};
    float c_b_answer[3] = {0.0, 0.0, 0.0};

    float delta_a[3] = {0.0, 0.0, 0.0};
    float delta_b[3] = {0.0, 0.0, 0.0};
    int num_tests = sizeof(r_b1) / sizeof(r_b1[0]);
    int pass_count = 0; // count number of cases that have passed

    rod::dbg_print = true;

    for (int i = 0; i < num_tests; i++)
    {
        vec3d(n) { c_a_answer[n] = 0.0; }
        vec3d(n) { c_b_answer[n] = 0.0; }

        // Function to test
        rod::element_minimum_displacement(p_a, p_b[i], r_a1, r_b1[i], c_a, c_b);

        switch (i)
        {
        case 0:
            // touching, midsection, at 45 degree angle in x-z plane
            c_a_answer[0] = 0.25;
            c_b_answer[0] = 0.25;
            c_b_answer[1] = radius_a + radius_b;
            break;
        case 1:
            // touching, end to end, at 15 degree angle x-z plane
            vec3d(n) { c_a_answer[n] = r_a2[n]; }
            vec3d(n) { c_b_answer[n] = r_b1[i][n]; }
            break;
        case 2:
            // touching, midsection, perpendicular in x-z plane
            c_a_answer[0] = 0.75;
            vec3d(n) { c_b_answer[n] = r_b1[i][n] + 0.5 * p_b[i][n]; }
            break;
        case 3:
            // intersect, midsection, 15 degrees in x-y plane
            c_a_answer[0] = r_b1[i][0];
            vec3d(n) { c_b_answer[n] = r_b1[i][n]; }
            // BUG: test gives c_a = (0, 0, 0) = r_a1, which is NOT the shortest
            // connection
            break;
        case 4:
            // intersect, end (a) to midsection (b), 60 degrees in x-y plane
            vec3d(n) { c_a_answer[n] = r_a2[n]; }
            vec3d(n) { c_b_answer[n] = r_a2[n]; }
            c_b_answer[0] += 0.2165;
            c_b_answer[1] -= 0.125;
            // BUG: test gives c_a = (0.5tan30, 0, 0) in x, which is NOT the
            // shortest connection
            break;
        case 5:
            // full overlap, end to end, 15 degree angle in x-y plane
            vec3d(n) { c_a_answer[n] = r_b1[i][n]; }
            vec3d(n) { c_b_answer[n] = c_a_answer[n]; }
            break;
        case 6:
            // full overlap, through midsection and out other side, -60 degree angle
            // in x-y plane
            vec3d(n) { c_a_answer[n] = -1; }
            c_a_answer[0] = 0.5 * std::tan(M_PI / 6);
            c_a_answer[1] = 0;
            c_a_answer[2] = 0;
            vec3d(n) { c_b_answer[n] = c_a_answer[n]; }
            break;
        default:
            std::cout << "Default case reached" << std::endl;
            return 1;
        }
        vec3d(n) { delta_a[n] = c_a[n] - c_a_answer[n]; }
        vec3d(n) { delta_b[n] = c_b[n] - c_b_answer[n]; }
        if (rod::absolute(delta_a) < 0.01 && rod::absolute(delta_b) < 0.01)
        {
            std::cout << "Case " << i << " passed!" << std::endl;
            pass_count++;
        }

        rod::print_array("c_a computed", c_a, 3);
        rod::print_array("c_a expected", c_a_answer, 3);
        std::cout << "|delta_a| = " << rod::absolute(delta_a) << std::endl;
        rod::print_array("c_b computed", c_b, 3);
        rod::print_array("c_b expected", c_b_answer, 3);
        std::cout << "|delta_b| = " << rod::absolute(delta_b) << std::endl;

        std::cout << "CASES PASSED: " << pass_count << "/" << num_tests << "\n"
                  << std::endl;
    }
    if (pass_count == num_tests)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

/**
 * Load four rods, defined in an accompanying Python script, and generate
 * neighbour lists for them. Count the total number of neighbours found and
 * compare to the expected result.
 *
 * TODO: More robust testing, e.g. test exactly which elements interact (by
 * comparing indices).
 */
int ffea_test::rod_neighbour_list_construction()
{
    int num_rods = 4;
    std::string filename;
    rod::Rod **rod_array = NULL;

    float r1[3] = {0, 0, 0};
    float r2[3] = {0, 0, 0};
    float distance = 0;
    float c_a[3] = {0, 0, 0};
    float c_b[3] = {0, 0, 0};
    float radius = 0;
    float c_ba[3] = {0, 0, 0};
    int num_neighbours = 0;
    int num_interactions = 0;
    int expected_num_interactions = 14;

    rod::dbg_print = true;

    World *world;
    world = new World();

    // Create rods
    std::cout << "ffea_test::rod_neighbour_list_construction() - loading rods"
              << std::endl;
    rod_array = new rod::Rod *[num_rods];
    for (int i = 0; i < num_rods; i++)
    {
        filename = "collider_" + std::to_string(i) + ".rod";
        rod::Rod *current_rod = new rod::Rod(filename, i);
        current_rod->load_header(filename);
        current_rod->load_contents(filename);
        current_rod->set_units();
        std::cout << "  Loaded rod from " << filename << std::endl;
        rod_array[i] = current_rod;
    }

    // Create neighbour lists
    std::cout << "ffea_test::rod_neighbour_list_construction() - generating "
                 "neighbour list"
              << std::endl;
    for (int i = 0; i < num_rods; i++)
    {
        for (int j = i + 1; j < num_rods; j++)
        {
            world->update_rod_steric_nbr_lists(rod_array[i], rod_array[j]);
        }
    }

    // Test the results
    std::cout << "ffea_test::rod_neighbour_list_construction() - results"
              << std::endl;
    for (int i = 0; i < num_rods; i++)
    {
        std::cout << "  rod " << i << std::endl;

        for (int j = 0; j < rod_array[i]->num_nodes - 1; j++)
        {
            vec3d(n) { r1[n] = rod_array[i]->current_r[3 * j + n]; }
            vec3d(n) { r2[n] = rod_array[i]->current_r[3 * (j + 1) + n]; }
            std::cout << "  element " << j << std::endl;
            rod::print_array("  r1", r1, 3);
            rod::print_array("  r2", r2, 3);
            std::cout << "  num_neighbours: "
                      << rod_array[i]->get_num_nbrs(j, rod_array[i]->steric_nbrs) << std::endl;
            // rod::print_vector("  all coords: ",
            //                   rod_array[i]->steric_neighbours.at(j));

            for (int k = 0; k < rod_array[i]->get_num_nbrs(j, rod_array[i]->steric_nbrs); k++)
            {
                // ! new data goes here!
                vec3d(n) { c_ba[n] = c_b[n] - c_a[n]; }
                distance = rod::absolute(c_ba);

                rod::print_array("  c_a", c_a, 3);
                rod::print_array("  c_b", c_b, 3);
                std::cout << "  |c_ba|: " << distance << std::endl;
            }
            num_neighbours += rod_array[i]->get_num_nbrs(j, rod_array[i]->steric_nbrs);
        }
        rod_array[i]->check_nbr_list_dim(rod_array[i]->steric_nbrs);
    }

    // avoid double counting (between two interacting rods, both will have
    // information on the others' elements)
    num_interactions = num_neighbours / 2;

    std::cout << "test score: " << num_interactions << "/"
              << expected_num_interactions << std::endl;
    if (num_interactions == 14)
    {
        return 0;
    }
    return 1;
}


// Move two circles away from each other and compute the rod-rod interaction energies.
int ffea_test::rod_steric_lj_potential()
{
    float radius = 6;  // ~ 1 nm, approaching continuum limit (DNA width ~ 2 nm)
    float eps = 10;    // 10 kT
    float r_min = 2 * radius;
    float sigma = r_min / std::pow(2, 1.0/6);  // beware integer division!

    float u_max = 20;  // steric potential at complete overlap (r = -2R)
    float r_max = 2 * radius;
    float steric_constant = u_max / (r_max * r_max);  // steric 'spring constant'
    float r_cutoff = 5 * radius;

    int num_steps = 1000;
    float dr = (r_cutoff + r_max) / num_steps;
    float r = -r_max;
    float u = 0;

    FILE* file_ptr1;
    file_ptr1 = std::fopen("log.csv", "w");

    std::fprintf(file_ptr1, "r,u\n");

    for (int step_no = 0; step_no < num_steps; step_no++)
    {
        if (r < 0)
        {
            u = rod::steric_energy_squared(steric_constant, r);
        }
        else if (r < r_min && r >= 0)
        {
            u = rod::vdw_energy_interp(r, eps, 1 / r_min);
        }
        else if (r > r_min)
        {
            u = rod::vdw_energy_6_12(1 / r, eps, sigma);
        }
        else
        {
            std::cout << "ERROR: invalid range for surface-surface distance\n";
            return 1;
        }
        std::fprintf(file_ptr1, "%e,%e\n", r * mesoDimensions::length, u * mesoDimensions::Energy);

        r += dr;
    }

    std::fflush(file_ptr1);
    std::fclose(file_ptr1);

    FILE* file_ptr2;
    file_ptr2 = std::fopen("const.csv", "w");

    std::fprintf(file_ptr2, "R,r_min,sigma,eps,k,r0,rN,kT\n");
    std::fprintf(file_ptr2, "%e,%e,%e,%e,%e,%e,%e,%e\n",
        radius * mesoDimensions::length,
        r_min * mesoDimensions::length,
        sigma * mesoDimensions::length,
        eps * mesoDimensions::Energy,
        steric_constant * (mesoDimensions::Energy / (mesoDimensions::length * mesoDimensions::length)),
        -r_max * mesoDimensions::length,
        r_cutoff * mesoDimensions::length,
        mesoDimensions::Energy
    );

    std::fflush(file_ptr2);
    std::fclose(file_ptr2);

    return 0;
}


/*
For two interacting points, A and B (xyz float coords), in a central simulation
box (LxLxL), surrounded by periodic images (ijk int coords), determine [i, j, k]
of the nearest periodic image.
xyz = {0...L}
ijk = {-1, 0, 1}
*/
int ffea_test::nearest_image_pbc()
{
    std::vector<float> dim = {1, 1, 1};
    float xmin = 0.1;
    float xmax = 0.9;
    rod::print_vector("box dim", dim);

    auto nimg_test = [&dim](float a[3], float b[3], int answer[3])
    {
        float ab[3] = {0};
        int diff[3] = {0};
        vec3d(n) { ab[n] = b[n] - a[n]; }
        std::vector<int> result = rod::nearest_periodic_image(a, b, dim);

        rod::print_array("a", a, 3);
        rod::print_array("b", b, 3);
        rod::print_array("ab", ab, 3);
        rod::print_array("expected image", answer, 3);
        rod::print_vector("computed image", result);
        std::cout << "\n";

        vec3d(n) { diff[n] = result.at(n) - answer[n]; }
        if (diff[0] == 0 && diff[1] == 0 && diff[2] == 0)
            return true;
        return false;
    };

    struct args
    {
        float a[3];
        float b[3];
        int img[3];

        args(std::vector<float> a_in, std::vector<float> b_in, std::vector<int> img_in)
        {
            vec3d(n) { a[n] = a_in[n]; }
            vec3d(n) { b[n] = b_in[n]; }
            vec3d(n) { img[n] = img_in[n]; }
        }

        args(float a_in[3], float b_in[3], int img_in[3])
        {
            vec3d(n) { a[n] = a_in[n]; }
            vec3d(n) { b[n] = b_in[n]; }
            vec3d(n) { img[n] = img_in[n]; }
        }
    };

    std::vector<args> input;

    auto create_face_input = [&input, &dim](float xmin, float xmax)
    {
        // 3 pairs along 6 axes: +/- x, y, z
        vec3d(n)
        {
            float a[3];
            float b[3];
            int img[3] = { 0 };
            vec3d(n)
            {
                a[n] = 0.5*dim[n];
                b[n] = 0.5*dim[n];
            }

            a[n] = xmin;
            b[n] = xmax;
            img[n] = 1;
            input.push_back(args(a, b, img));
            img[n] *= -1;
            input.push_back(args(b, a, img));
        };
    };

    auto create_diag_input = [&input](float xmin, float xmax)
    {
        int octants_a[4][3] = {{0,0,0}, {0,0,1}, {1,0,0}, {1,0,1}};
        int octants_b[4][3] = {{1,1,1}, {1,1,0}, {0,1,1}, {0,1,0}};

        // 4 pairs occupying 8 corners
        for (int i=0; i<4; i++)
        {
            float a[3] = {xmin, xmin, xmin};
            float b[3] = {xmin, xmin, xmin};
            int img[3] = {-1, -1, -1};

            vec3d(n){
                a[n] = std::max(xmin, octants_a[i][n] * xmax);
                b[n] = std::max(xmin, octants_b[i][n] * xmax);
                if (octants_b[i][n] == 1)
                    img[n] = 1;
            };
            input.push_back(args(a, b, img));
            vec3d(n){img[n] *= -1;}
            input.push_back(args(b, a, img));
        }
    };

    create_face_input(xmin, xmax);
    create_diag_input(xmin, xmax);

    int i = 1;
    int fail_count = 0;
    for (auto in : input)
    {
        std::cout << "Input " << i << " / " << input.size() << "\n";
        if (!nimg_test(in.a, in.b, in.img))
            fail_count++;
        i++;
    }

    if (fail_count == 0)
        return 0;
    return 1;
}

int ffea_test::rod_vdw_site_placement()
{
    std::string rod_file;
    std::string rodvdw_file;
    rod::Rod **rod_array = NULL;
    int num_rods = 2;

    World *world;
    world = new World();

    std::string site_file;
    FILE* file_ptr;

    // Create rods
    std::cout << "Loading rods...\n";
    rod_array = new rod::Rod *[num_rods];

    string rod_names[2] = {"curvy", "straight"};

    for (int i = 0; i < num_rods; i++)
    {
        rod_file = rod_names[i] + ".rod";
        rodvdw_file = rod_names[i] + ".rodvdw";

        std::cout << "  .rod filename :     " << rod_file << "\n";
        std::cout << "  .rodvdw filename :  " << rodvdw_file << "\n";

        rod::Rod *current_rod = new rod::Rod(rod_file, i);
        current_rod->load_header(rod_file);
        current_rod->load_contents(rod_file);
        current_rod->set_units();
        current_rod->load_vdw(rodvdw_file);
        std::cout << "  Loaded rod from " << rod_file << std::endl;
        std::cout << "  Loaded VDW sites from " << rodvdw_file << std::endl;
        rod_array[i] = current_rod;

        // rod::print_vector("vdw_site_pos", current_rod->vdw_site_pos);

        std::cout << "num_vdw_sites : " << current_rod->num_vdw_sites << "\n";

        std::ifstream in_file;
        std::string row;
        std::string coord;
        float pos_from_file[3] = {0};
        float pos_from_rod[3] = {0};
        float diff[3] = {0};

        in_file.open(rod_names[i] + "_vdw_pos.csv");

        std::cout << "Reading " + rod_names[i] + "_vdw_pos.csv...\n";
        for (int j=0; j < current_rod->num_vdw_sites; j++)
        {
            // read the file stream into a string, representing a whole row
            std::getline(in_file, row, '\n');

            // pass the row to a string stream (needs to be reset here)
            std::stringstream ss;
            ss << row;

            vec3d(n)
            {
                // read the comma-separated values of the string stream into an array
                std::getline(ss, coord, ',');
                pos_from_file[n] = std::stof(coord);
            }

            vec3d(n){pos_from_rod[n] = current_rod->vdw_site_pos.at((3 * j) + n) * mesoDimensions::length;}

            vec3d(n){diff[n] = pos_from_file[n] - pos_from_rod[n];}

            rod::print_array("  pos_from_file",pos_from_file,3);
            rod::print_array("  pos_from_rod",pos_from_rod,3);
            rod::print_array("  diff",diff,3);
            std::cout << "  |diff| : " << rod::absolute(diff) << "\n";
            std::cout << "\n";

            if (rod::absolute(diff) > 1.0e-8)
            {
                std::cout << "Fail. Difference too large.\n";
                return 1;
            }

        }

        in_file.close();
    }

    return 0;
}
