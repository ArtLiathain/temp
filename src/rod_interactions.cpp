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
 *      rod_vdw.cpp
 *	Author: Ryan Cocking, University of Leeds
 *	Email: bsrctb@leeds.ac.uk
 */

#include "rod_interactions.h"

namespace rod
{

/*
================================================================================
    INTERACTION OBJECT
================================================================================
*/

// Steric interactions
InteractionData::InteractionData(int rod_id_a, int rod_id_b, int elem_id_a,
    int elem_id_b, float radius_a, float radius_b, float c_a[3], float c_b[3],
    float shift[3], float r_a[3], float r_b[3])
{
    rod_id_self = rod_id_a;
    rod_id_nbr = rod_id_b;
    elem_id_self = elem_id_a;
    elem_id_nbr = elem_id_b;
    radius_self = radius_a;
    radius_nbr = radius_b;
    vec3d(n){contact_self[n] = c_a[n];}
    vec3d(n){contact_nbr[n] = c_b[n];}
    vec3d(n){c_ab[n] = c_b[n] - c_a[n];}
    vec3d(n){img_shift[n] = shift[n];}
    vec3d(n){r_self[n] = r_a[n];}
    vec3d(n){r_nbr[n] = r_b[n];}

    epsilon = 0;
    sigma = 0;
    r_min = 0;
    r_min_inv = 0;
}

// VDW interactions
InteractionData::InteractionData(int rod_id_a, int rod_id_b, int elem_id_a,
    int elem_id_b, float radius_a, float radius_b, float c_a[3], float c_b[3],
    float shift[3], float r_a[3], float r_b[3], float eps, float sig)
{
    rod_id_self = rod_id_a;
    rod_id_nbr = rod_id_b;
    elem_id_self = elem_id_a;
    elem_id_nbr = elem_id_b;
    radius_self = radius_a;
    radius_nbr = radius_b;
    vec3d(n){contact_self[n] = c_a[n];}
    vec3d(n){contact_nbr[n] = c_b[n];}
    vec3d(n){c_ab[n] = c_b[n] - c_a[n];}
    vec3d(n){img_shift[n] = shift[n];}
    vec3d(n){r_self[n] = r_a[n];}
    vec3d(n){r_nbr[n] = r_b[n];}
    epsilon = eps;
    sigma = sig;
    r_min = std::pow(2, 1./6) * sigma;
    r_min_inv = 1 / r_min;

}

bool InteractionData::elements_intersect()
{
    float c_ab[3] = {
        this->contact_nbr[0] - this->contact_self[0],
        this->contact_nbr[1] - this->contact_self[1],
        this->contact_nbr[2] - this->contact_self[2]
    };

    if (rod::dbg_print)
    {
        std::cout << "|c_ab|: " << rod::absolute(c_ab) << "\n";
        std::cout << "radius sum: " << this->radius_self + this->radius_nbr << "\n";
    }

    if (rod::absolute(c_ab) < this->radius_self + this->radius_nbr)
        return true;
    else
        return false;
}

/*
================================================================================
    STERIC INTERACTIONS
================================================================================
*/

/**
 * If an interaction point (c) lies outside the finite length of its rod element,
 * replace its value with the nearest node. Required for certain situations, e.g.
 * nearly-parallel rods.
 */
void finite_length_correction(float c[3], float r[3], float p[3])
{
    float r_c[3] = { 0 };
    float dot = 0;

    vec3d(n) { r_c[n] = c[n] - r[n]; }
    dot = dot_product_3x1(p, r_c);

    if (dot <= 0)
    {
        vec3d(n) { c[n] = r[n]; }
    }
    else if (dot >= rod::absolute(p) * rod::absolute(p))
    {
        vec3d(n) { c[n] = r[n] + p[n]; }
    }

    // if (rod::dbg_print)
    // {
    //     std::cout << "finite_length_correction():\n";
    //     printf("  p.(c - r) : %.3e\n", dot);
    // }
}

// stupid implementation to get around a linker error
void finite_length_correction(float c[3], float r[3], float p[3], OUT float c_out[3])
{
    float r_c[3] = { 0 };
    float dot = 0;

    vec3d(n) { r_c[n] = c[n] - r[n]; }
    dot = dot_product_3x1(p, r_c);

    if (dot <= 0)
    {
        vec3d(n) { c_out[n] = r[n]; }
    }
    else if (dot >= rod::absolute(p) * rod::absolute(p))
    {
        vec3d(n) { c_out[n] = r[n] + p[n]; }
    }

    // if (rod::dbg_print)
    // {
    //     std::cout << "finite_length_correction():\n";
    //     printf("  p.(c - r) : %.3e\n", dot);
    // }
}


/** Compare the centreline displacement, c_ab, to the following
 * rod-rod displacements:
 *
 * d1 = c_b - r_a1
 * d2 = c_b - r_a2
 * d3 = c_a - r_b1
 * d4 = c_a - r_b2
 *
 * Find the smallest displacement and, if needed, adjust c_ab such that it
 * becomes the smallest displacement.
*/

void nearest_node_correction(float c_a[3], float c_b[3],
    float r_a[3], float r_b[3], float p_a[3], float p_b[3])
{
    float displ[3][3] = { 0 };
    float r_a2[3] = { 0 };
    float r_b2[3] = { 0 };
    std::vector<float> dist(3, 0);

    vec3d(n) { r_a2[n] = r_a[n] + p_a[n]; }
    vec3d(n) { r_b2[n] = r_b[n] + p_b[n]; }
    vec3d(n) { displ[0][n] = c_b[n] - c_a[n]; }
    vec3d(n) { displ[1][n] = c_b[n] - r_a[n]; }
    vec3d(n) { displ[2][n] = c_b[n] - r_a2[n]; }

    // index corresponding to the smallest distance
    vec3d(n) { dist[n] = rod::absolute(displ[n]); }
    auto iter = std::min_element(dist.begin(), dist.end());
    int index_min = std::distance(dist.begin(), iter);

    if (index_min == 1)
        vec3d(n) { c_a[n] = r_a[n]; }
    else if (index_min == 2)
        vec3d(n) { c_a[n] = r_a2[n]; }
    else if (index_min != 0)
        throw std::out_of_range("Invalid min index in rod::nearest_node_connection()");

    // if (rod::dbg_print)
    // {
    //     std::cout << "nearest_node_correction():\n";
    //     print_array("  r_a1", r_a, 3);
    //     print_array("  r_a2", r_a2, 3);
    //     printf("  |c_ab| :       %.3f\n", dist[0]);
    //     printf("  |c_b - r_a1| : %.3f\n", dist[1]);
    //     printf("  |c_b - r_a2| : %.3f\n", dist[2]);
    //     print_array("  c_a (corrected)", c_a, 3);
    //     std::cout << "\n";
    // }
}

/** @brief Compute the two points, c_a and c_b, that form the centreline
 * displacement (interaction vector) joining two rod elements together, where
 * c_a sits on the element p_a.
 *
 * The displacement for oblique rod elements is calculated using an infinite
 * line assumption, so a correction is applied to account for the finite length
 * of elements.
*/
void element_minimum_displacement(float p_a[3], float p_b[3], float r_a[3],
    float r_b[3], float c_a[3], float c_b[3])
{
    float l_a[3] = {0}; // l_a = p_a / |p_a|
    float l_b[3] = {0};
    float l_a_cross_l_b[3] = {0};
    float n_a[3] = {0};
    float n_b[3] = {0};
    float r_ab[3] = {0};
    float r_ba[3] = {0};

    normalize(p_a, l_a);
    normalize(p_b, l_b);
    cross_product(l_a, l_b, l_a_cross_l_b);

    // oblique / perpendicular elements
    if (rod::absolute(l_a_cross_l_b) > 1e-5)
    {
        cross_product(l_a, l_a_cross_l_b, n_a);
        cross_product(l_b, l_a_cross_l_b, n_b);

        vec3d(n) { r_ab[n] = r_b[n] - r_a[n]; }
        vec3d(n) { r_ba[n] = -r_ab[n]; }

        vec3d(n)
        {
            c_a[n] = r_a[n] +
                dot_product_3x1(r_ab, n_b) / dot_product_3x1(l_a, n_b) * l_a[n];
        }
        vec3d(n)
        {
            c_b[n] = r_b[n] +
                dot_product_3x1(r_ba, n_a) / dot_product_3x1(l_b, n_a) * l_b[n];
        }
    }
    // parallel elements
    else
    {
        vec3d(n) { c_a[n] = r_a[n] + 0.5 * p_a[n]; }
        vec3d(n) { c_b[n] = r_b[n] + 0.5 * p_b[n]; }
    }

    // if (rod::dbg_print)
    // {
    //     std::cout << "minimum distance between rod elements:\n";
    //     print_array("  p_a", p_a, 3);
    //     print_array("  p_b", p_b, 3);
    //     print_array("  l_a", l_a, 3);
    //     print_array("  l_b", l_b, 3);
    //     print_array("  l_a x l_b", l_a_cross_l_b, 3);
    //     std::cout << "  |l_a x l_b|: " << rod::absolute(l_a_cross_l_b) << "\n";
    //     print_array("  n_b", n_a, 3);
    //     print_array("  n_b", n_b, 3);
    //     print_array("  r_a", r_a, 3);
    //     print_array("  r_b", r_b, 3);
    //     print_array("  r_ab", r_ab, 3);
    //     print_array("  r_ba", r_ba, 3);
    //     print_array("  c_a (initial)", c_a, 3);
    //     print_array("  c_b (initial)", c_b, 3);
    //     std::cout << "\n";
    // }

    // enforce finite length of rod element
    finite_length_correction(c_a, r_a, p_a);
    finite_length_correction(c_b, r_b, p_b);

    // if a node is closer than c_a or c_b, use that instead
    nearest_node_correction(c_a, c_b, r_a, r_b, p_a, p_b);
    nearest_node_correction(c_b, c_a, r_b, r_a, p_b, p_a);
}

/*
Return the integer coordinates of the nearest image with which to correct the
displacement of an interacting element pair.
For displacment pointing from A to B, the returned image is that nearest to B,
so the corrected displacement would require subtracting a box length from B
*/
std::vector<int> nearest_periodic_image(float a[3], float b[3], std::vector<float> box_dim)
{
    std::vector<int> img(3, 0);
    vec3d(n) { img.at(n) = std::floor((b[n] - a[n] + 0.5 * box_dim[n]) / box_dim[n]); }

    if (dbg_print)
    {
        std::cout << "nearest_periodic_image:\n";
        print_array("  a", a, 3);
        print_array("  b", b, 3);
        print_vector("  box_dim", box_dim);
        print_vector("  img", img);
    }

    return img;
}

/*
Check if two rod elements, a and b, interact by calculating the shortest
distance between them and comparing to the sum of their radii. If this
passes, the interaction information is added to both elements' neighbour
lists.
*/
void set_steric_nbrs(int rod_id_a, int rod_id_b, int elem_id_a,
    int elem_id_b, float p_a[3], float p_b[3], float r_a[3], float r_b[3],
    float radius_a, float radius_b, std::vector<InteractionData> &neighbours_a,
    std::vector<InteractionData> &neighbours_b, bool periodic, std::vector<float> box_dim)
{
    float mid_a[3] = {0};
    float mid_b[3] = {0};
    float mid_ab[3] = {0};
    float c_a[3] = {0};
    float c_b[3] = {0};
    float c_ab[3] = {0};
    std::vector<int> img = {0, 0, 0};
    float shift[3] = {0};
    float radius_sum = radius_a + radius_b;

    rod::get_element_midpoint(p_a, r_a, mid_a);
    rod::get_element_midpoint(p_b, r_b, mid_b);
    vec3d(n){mid_ab[n] = mid_b[n] - mid_a[n];}

    // PBC correction
    if (periodic)
    {
        img = nearest_periodic_image(mid_a, mid_b, box_dim);
        vec3d(n){shift[n] = box_dim[n] * img.at(n);}
        vec3d(n){r_b[n] -= shift[n];}
    }

    rod::get_element_midpoint(p_b, r_b, mid_b);
    vec3d(n){mid_ab[n] = mid_b[n] - mid_a[n];}

    if (rod::dbg_print)
    {
        std::cout << "cull distant elements:\n";
        std::printf("  elem index a :  %d\n", elem_id_a);
        std::printf("  elem index b :  %d\n", elem_id_b);
        std::printf("  |p_a| :         %.3f\n", rod::absolute(p_a));
        std::printf("  |p_b| :         %.3f\n", rod::absolute(p_b));
        std::printf("  radius sum :    %.3f\n", radius_sum);
        std::printf("  |midpoint ab| : %.3f\n", rod::absolute(mid_ab));
        std::printf("  cutoff     :    %.3f\n", std::max(rod::absolute(p_a), rod::absolute(p_b)) + radius_sum);
    }

    float cutoff = std::max(rod::absolute(p_a), rod::absolute(p_b)) + radius_sum;

    if (rod::absolute(mid_ab) < cutoff)
    {
        rod::element_minimum_displacement(p_a, p_b, r_a, r_b, c_a, c_b);
        vec3d(n) { c_ab[n] = c_b[n] - c_a[n]; }

        if (rod::dbg_print)
        {
            std::cout << "rod-rod distance:\n";
            printf("  |c_ab| : %.3e\n", rod::absolute(c_ab));
        }

        if (rod::absolute(c_ab) > 1e-5 && rod::absolute(c_ab) < radius_sum)
        {
            if (rod::dbg_print)
            {
                std::cout << "  neighbour pair: " << rod_id_a << "|" << elem_id_a
                << " and " << rod_id_b << "|" << elem_id_b << "\n\n";
                std::cout << "  generating interaction structs\n";
            }

            InteractionData stericDataA(
                rod_id_a,
                rod_id_b,
                elem_id_a,
                elem_id_b,
                radius_a,
                radius_b,
                c_a,
                c_b,
                shift,
                r_a,
                r_b);
            neighbours_a.push_back(stericDataA);

            InteractionData stericDataB(
                rod_id_b,
                rod_id_a,
                elem_id_b,
                elem_id_a,
                radius_b,
                radius_a,
                c_b,
                c_a,
                shift,
                r_b,
                r_a);
            neighbours_b.push_back(stericDataB);
        }
        else if (rod::absolute(c_ab) <= 1e-5)
        {
            std::cout << "ERROR!\n";
            throw std::runtime_error("Rod elements have fully overlapped.");
        }
    }
    else if (rod::dbg_print && rod::absolute(mid_ab) >= cutoff)
    {
        std::cout << "  ignored; outside steric regime\n\n";
    }
    else if (rod::dbg_print)
    {
        throw std::runtime_error("Invalid distance to rod::set_steric_nbrs()");
    }
}

float steric_energy_squared(float force_scaling_factor, float intersect_distance)
{
    return force_scaling_factor * intersect_distance * intersect_distance;
}

// Perturbation in +c_ab/-c_ab, where c_ab points away from current element
float intersection_distance(float delta, float c_a[3], float c_b[3], float radius_sum)
{
    float c_ab[3] = { 0 };
    float c_ab_norm[3] = { 0 };
    float c_ab_ptb[3] = { 0 };

    vec3d(n) { c_ab[n] = c_b[n] - c_a[n]; }
    rod::normalize(c_ab, c_ab_norm);
    vec3d(n) { c_ab_ptb[n] = c_ab[n] + delta * c_ab_norm[n]; }

    return rod::absolute(c_ab_ptb) - radius_sum;
}

// No perturbation
float intersection_distance(float c_a[3], float c_b[3], float radius_sum)
{
    float c_ab[3] = { 0 };
    vec3d(n) { c_ab[n] = c_a[n] - c_b[n]; }

    return rod::absolute(c_ab) - radius_sum;
}

// Return the steric collision energy of a rod element as a vector with the
// format [r+dr, r-dr, r]. Perturbation is applied along the intersection vector.
std::vector<float> element_steric_energy(float delta, float force_constant,
    float radius_sum, float c_a[3], float c_b[3])
{
    float dist[3] = {0};
    float energy[3] = {0};

    dist[0] = intersection_distance(0.5 * delta, c_a, c_b, radius_sum);
    dist[1] = intersection_distance(-0.5 * delta, c_a, c_b, radius_sum);
    dist[2] = intersection_distance(c_a, c_b, radius_sum);

    energy[0] = steric_energy_squared(force_constant, dist[0]);
    energy[1] = steric_energy_squared(force_constant, dist[1]);
    energy[2] = steric_energy_squared(force_constant, dist[2]);

    if (rod::dbg_print)
    {
        float c_ab[3] = { 0 };
        vec3d(n) { c_ab[n] = c_b[n] - c_a[n]; }
        // std::cout << "element steric energy:\n";
        // std::cout << "  delta : " << delta << "\n";
        // std::cout << "  force_constant : " << force_constant << "\n";
        // std::cout << "  radius_sum : " << radius_sum << "\n";
        // print_array("  c_a", c_a, 3);
        // print_array("  c_b", c_b, 3);
        // std::cout << "  |c_ab| : " << rod::absolute(c_ab) << "\n";
        // printf("  perturbed distance +ve : %.3e\n", dist[0]);
        // printf("  perturbed distance -ve : %.3e\n", dist[1]);
        // printf("  distance : %.3e\n", dist[2]);
        // printf("  perturbed energy +ve : %.3e\n", energy[0]);
        // printf("  perturbed energy -ve : %.3e\n", energy[1]);
        // printf("  energy : %.3e\n", energy[2]);
        // std::cout << "\n";

        if (rod::absolute(c_ab) > radius_sum)
        {
            std::cout << "ERROR!\n";
            throw std::runtime_error("Centreline distance cannot be "
                "greater than radius sum when calculating steric energy.");
        }
    }

    return std::vector<float> {energy[0], energy[1], energy[2]};
}

// Interpolate the force applied to an element onto both nodes
std::vector<float> node_force_interpolation(float contact[3], float node_start[3],
    float element_length, const std::vector<float> &element_force)
{
    float displacement[3] = { 0 };
    float l1 = 0;
    float l2 = 0;
    std::vector<float> force(6, 0);  // x0, y0, z0, x1, y1, z1

    vec3d(n) { displacement[n] = contact[n] - node_start[n]; }

    l1 = rod::absolute(displacement);
    l2 = element_length - l1;

    vec3d(n) { force[n] = element_force[n] * (element_length - l1) / element_length; }
    vec3d(n) { force[n + 3] = element_force[n] * (element_length - l2) / element_length; }

    if (rod::dbg_print)
    {
        std::cout << "node force interpolation:\n";
        print_array("  contact", contact, 3);
        print_array("  node_start", node_start, 3);
        std::cout << "  element_length, L : " << element_length << "\n";
        print_vector("  element_force", element_force);
        print_array("  displacement", displacement, 3);
        std::cout << "  l1 : " << l1 << "\n";
        std::cout << "  l2 : " << l2 << "\n";
        std::cout << "  weight start node : " << (element_length - l1) / element_length << "\n";
        std::cout << "  weight end node : " << (element_length - l2) / element_length << "\n";
        print_vector("  force start node", std::vector<float>(force.begin(), std::next(force.begin(), 3)) );
        print_vector("  force end node", std::vector<float>(std::next(force.begin(), 3), force.end()) );
        std::cout << "\n";
    }

    if (l1 > element_length)
    {
        std::string err = "Distance along rod element (" + std::to_string(l1) +
            ") cannot be greater than the element length (" +
            std::to_string(element_length) + ")\n";
        throw std::invalid_argument(err);
    }

    return force;
}

// ! - should ideally be in their own file
/*
================================================================================
    VAN DER WAALS INTERACTIONS
================================================================================
*/

float vdw_energy_6_12(float r_inv, float eps, float sig)
{
    float sr = sig * r_inv;
    float sr3 = sr * sr * sr;
    float sr6 = sr3 * sr3;
    float sr12 = sr6 * sr6;
    return 4 * eps * (sr12 - sr6);
}

float vdw_force_6_12(float r_inv, float eps, float sig)
{
    float sr = sig * r_inv;
    float sr3 = sr * sr * sr;
    float sr6 = sr3 * sr3;
    float sr12 = sr6 * sr6;
    return -24 * eps * r_inv * (2 * sr12 - sr6);
}

// r_min_inv = 1 / (2^1/6 * sigma)
// r = surface-surface distance
float vdw_energy_interp(float r, float eps, float r_min_inv)
{
    float rr = r * r_min_inv;
    return eps * (2* rr * rr * rr - 3 * rr * rr);
}

// r_min_inv = 1 / (2^1/6 * sigma)
float vdw_force_interp(float r, float eps, float r_min_inv)
{
    float rr = r * r_min_inv;
    return 6 * eps * r_min_inv * (rr * rr - rr);
}

VDWSite::VDWSite(const int rodid, const int siteid, const int vdwtype, const float lrod,
    const float *r_rod, const float contour_length, const int num_nodes)
{
    rod_id = rodid;
    site_id = siteid;
    vdw_type = vdwtype;
    L_rod = lrod;

    this->get_parent_element(r_rod, contour_length, L_rod, num_nodes);
    this->update_position(r_rod);
}

void VDWSite::print_info() const
{
    std::printf("VDWSite info:\n");
    std::printf("\trod_id:   %d\n", this->rod_id);
    std::printf("\telem_id:  %d\n", this->elem_id);
    std::printf("\tsite_id:  %d\n", this->site_id);
    std::printf("\tvdw_type: %d\n", this->vdw_type);
    std::printf("\tL_rod:    %.3f\n", this->L_rod);
    std::printf("\tL_elem:   %.3f\n", this->L_elem);
    rod::print_array("\tposition", this->pos, 3);
}

// Ensure parent element is set before this is called
void VDWSite::update_position(const float* r_rod)
{
    float r1[3] = {0};
    float r2[3] = {0};
    float p[3] = {0};
    vec3d(n)
    {
        r1[n] = r_rod[this->elem_id * 3 + n];
        r2[n] = r_rod[(this->elem_id + 1) * 3 + n];
    }
    rod::get_p_i(r1, r2, p);
    vec3d(n) { this->pos[n] = r1[n] + p[n] * this->L_elem; }

}


// Determine the rod element that a VDW site belongs to.
void VDWSite::get_parent_element(const float* r_rod, const float contour_length,
    const float norm_length_along_rod, const int num_nodes)
{
    float r1[3] = {0};
    float r2[3] = {0};
    float p[3] = {0};
    float psum = 0;
    float psum_prev = 0;
    float length_along_rod = norm_length_along_rod * contour_length;

    // traverse along the rod
    for (int node = 0; node < num_nodes - 1; ++node)
    {
        vec3d(n)
        {
            r1[n] = r_rod[node * 3 + n];
            r2[n] = r_rod[(node + 1) * 3 + n];
        }
        rod::get_p_i(r1, r2, p);

        float pmag = rod::absolute(p);
        psum += pmag;
        float norm_length_along_elem = (length_along_rod - psum_prev) / pmag;
        psum_prev = psum;

        // return if we overstep
        if (psum > length_along_rod)
        {
            this->elem_id = node;
            this->L_elem = norm_length_along_elem;
            return;
        }
    }

    // otherwise, set to final element
    this->elem_id = num_nodes - 2;
    this->L_elem = 1;
}


void set_vdw_nbrs(const VDWSite site_a, const VDWSite site_b, float p_a[3], float p_b[3],
    float r_a[3], float r_b[3], float radius_a, float radius_b,
    std::vector<InteractionData>& nbr_a, std::vector<InteractionData>& nbr_b,
    bool periodic, std::vector<float> box_dim, float vdw_cutoff,
    float epsilon, float sigma)
{
    float c_a[3] = { 0 };
    float c_b[3] = { 0 };
    float c_ab[3] = { 0 };
    std::vector<int> img = { 0, 0, 0 };
    float shift[3] = { 0 };
    float mag = 0;
    float radius_sum = radius_a + radius_b;

    // centreline distance is just point-to-point distance of sites, nothing fancy
    vec3d(n){c_a[n] = site_a.pos[n];}
    vec3d(n){c_b[n] = site_b.pos[n];}

    if (periodic)
    {
        img = nearest_periodic_image(c_a, c_b, box_dim);
        vec3d(n) { shift[n] = box_dim[n] * img.at(n); }
        vec3d(n) { c_b[n] -= shift[n]; }
    }

    vec3d(n) { c_ab[n] = c_b[n] - c_a[n]; }
    mag = rod::absolute(c_ab);

    if (rod::dbg_print)
    {
        std::printf("  |c_ab|             : %.3e\n", mag);
        std::printf("  surf-surf distance : %.3e\n", mag - radius_sum);
        std::printf("  vdw_cutoff         : %.3e\n", vdw_cutoff);
        std::printf("  eps                : %.3e\n", epsilon);
        std::printf("  sig                : %.3e\n", sigma);
        site_a.print_info();
        site_b.print_info();
    }

    // Interactions are calculated based on surface-surface distance
    if (mag - radius_sum > 0 && mag - radius_sum < vdw_cutoff)
    {
        InteractionData vdwDataA(
            site_a.rod_id,
            site_b.rod_id,
            site_a.elem_id,
            site_b.elem_id,
            radius_a,
            radius_b,
            c_a,
            c_b,
            shift,
            r_a,
            r_b,
            epsilon,
            sigma);
        nbr_a.push_back(vdwDataA);

        InteractionData vdwDataB(
            site_b.rod_id,
            site_a.rod_id,
            site_b.elem_id,
            site_a.elem_id,
            radius_b,
            radius_a,
            c_b,
            c_a,
            shift,
            r_b,
            r_a,
            epsilon,
            sigma);
        nbr_b.push_back(vdwDataB);
    }
    else if (rod::dbg_print && mag - radius_sum >= vdw_cutoff)
    {
        std::cout << "  ignored; exceeded vdw cutoff\n";
    }
    else if (rod::dbg_print && mag - radius_sum <= 0)
    {
        std::cout << "  ignored; within steric regime\n";
    }
    else if (rod::dbg_print)
    {
        throw std::runtime_error("Invalid distance to rod::set_vdw_nbrs()");
    }
}

//    __      _
//  o'')}____//  I AM DEBUG DOG. PUT ME IN YOUR
//   `_/      )  SOURCE CODE AND I WILL EAT THE
//   (_(_/-(_/   BUGS. WOOF WOOF!

} // namespace rod
