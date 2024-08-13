#include "rod_structure.h"
#include <string>
#include <cmath>

// A less weird way to access the contents of our arrays representing our vectors
static const int x = 0;
static const int y = 1;
static const int z = 2;

// For clarity, anytime you see an index [x], we're referring to the x
// dimension.

// Computing the energy for a perturbation requires four segments, defined
// by 5 nodes. They are i-2 to i+1, listed here.
static const int im2 = 0; ///< index of i-2nd thing
static const int im1 = 1; ///< index of i-1st thing
static const int i = 2; ///< index of ith thing
static const int ip1 = 3; ///< index of i+1th thing

float twist_energy_perturbation_test(
        rod::Rod rod_to_test,
        int perturbation_dimension,
        int p_i_node_no,
        float perturbation_amount
    ){
        
    //float perturbation_amount = rod_to_test.perturbation_amount;
    //int perturbation_dimension = rod_to_test.perturbation_dimension;
    //int start_cutoff = rod_to_test.start_cutoff;
    //int end_cutoff = rod_to_test.end_cutoff;
    //int p_i_node_node_no = rod_to_test.p_i_node_no;
        
    float *B_matrix = rod_to_test.B_matrix;
    float *material_params = rod_to_test.material_params;
    float *r_all = rod_to_test.current_r;
    float *r_all_equil = rod_to_test.equil_r;
    float *m_all = rod_to_test.current_m;
    float *m_all_equil = rod_to_test.equil_m;
    
    int start_cutoff;
    int end_cutoff;
    int *start_cutoff_ptr = &start_cutoff;
    int *end_cutoff_ptr = &end_cutoff; // for the multiple return values
    rod::set_cutoff_values(p_i_node_no, rod_to_test.num_nodes, start_cutoff_ptr, end_cutoff_ptr);
                
    //feenableexcept( FE_INVALID | FE_OVERFLOW );
        
    // Put a 5-node segment onto the stack.
    // We need to make a copy of it, because we'l be modifying it for our
    // Numerical differentiation later on.

    float B_equil[4][4];
    rod::load_B_all(B_equil, B_matrix, p_i_node_no);
    
    // We'll end up modifying this, but we need the original later to update the material frame
    float original_p[4][3];   
    rod::load_p(original_p, r_all, p_i_node_no);

    float p[4][3]; // the perturbed e
    float p_equil[4][3];
    float m[4][3];
    float m_equil[4][3];
    float material[4][3]; // 0 = k (stretch), 1 = beta (twist), 2 = unused (for now)
    
//    rod::print_array("rod_to_test.current_m", rod_to_test.current_m, 18);
    
    // Compute e from m, and load into an appropriate data structure (a 2-D array)
    rod::load_p(p, r_all, p_i_node_no);
    rod::load_p(p_equil, r_all_equil, p_i_node_no);
    rod::load_m(m, m_all, p_i_node_no);    
    rod::load_m(m_equil, m_all_equil, p_i_node_no);
    rod::load_m(material, material_params, p_i_node_no);
    
//    std::cout << "p_i = [";
//    for (int inc = 0; inc < 4; inc++){
//        rod::print_array("", p[inc], 3);
//    }
//    std::cout << "\n";
    
    // Normalize m1, just to be sure (the maths doesn't work if it's not normalized)
    rod::normalize_all(m);
    rod::normalize_all(m_equil);
    
//    std::cout << "p_i_norm = [";
//    for (int inc = 0; inc < 4; inc++){
//        rod::print_array("", p[inc], 3);
//    }
//    std::cout << "\n";
    
    std::cout << "orig p = " << p[i][x] << ", " << p[i][y] << ", " << p[i][z] << "\n";
    
    // Apply our perturbation in x, y or z (for numerical differentiation)
    if (perturbation_dimension < 4 && perturbation_amount != 0){ //e.g. if we're perturbing x, y, or z
        p[im1][perturbation_dimension] += perturbation_amount;
        p[i][perturbation_dimension] -= perturbation_amount;
    }
    
    // If we perturb our angle instead, we apply a rodrigues rotation.
    if(perturbation_dimension == 4 && perturbation_amount != 0){ // if we're perturbing the twist
        rod::rodrigues_rotation(m[i], p[i], perturbation_amount, m[i]);
    }
    
    std::cout << "orig m = " << m[i][x] << ", " << m[i][y] << ", " << m[i][z] << "\n";
    std::cout << " new p = " << p[i][x] << ", " << p[i][y] << ", " << p[i][z] << "\n";
    
    //std::cout << "     m = [";
    //for (int inc = 0; inc < 4; inc++){
    //    std::cout << "[ " << m[inc][x] << ", " << m[inc][y] << ", " << m[inc][z] << "] ";
    //}
    //std::cout << "]\n";
    //
    //If we've perturbed it in x, y, or z, we need to update m, and then adjust it to make sure it's perpendicular
    if (perturbation_dimension < 4 && perturbation_amount != 0){ 
        rod::update_m1_matrix_all(m, original_p, p, m, start_cutoff, end_cutoff);
    }
    
    std::cout << " new m = " << m[i][x] << ", " << m[i][y] << ", " << m[i][z] << "\n";
    
    //
    //std::cout << " new m = [";
    //for (int inc = 0; inc < 4; inc++){
    //    std::cout << "[ " << m[inc][x] << ", " << m[inc][y] << ", " << m[inc][z] << "] ";
    //}
    //std::cout << "]\n";
    //
    //std::cout << "     p = [";
    //for (int inc = 0; inc < 4; inc++){
    //    std::cout << "[ " << p[inc][x] << ", " << p[inc][y] << ", " << p[inc][z] << "] ";
    //}
    //std::cout << "]\n";
    
    //std::cout << "orig p = [";
    //for (int inc = 0; inc < 4; inc++){
    //    std::cout << "[ " << original_p[inc][x] << ", " << original_p[inc][y] << ", " << original_p[inc][z] << "] ";
    //}
    //std::cout << "]\n";
    
    //std::cout << "start_cutoff = " << start_cutoff << "\n";
    
    //float m_orig_temp[4][3];
    //rod::load_m(m_orig_temp, m_all, p_i_node_no);    
    //float m1_updated_twisten = rod::get_twist_energy(material[i][1], m[i], m[im1], m_orig_temp[i], m_orig_temp[im1], p[im1], p[i], original_p[im1], original_p[i]);
    //std::cout << "TWIST ENERGY DELIVERED BY UPDATING M1 (should be 0!): " << m1_updated_twisten << "\n";
    
//    std::cout << "p_i_new = [";
//    for (int inc = 0; inc < 4; inc++){
//        rod::print_array("", p[inc], 3);
//    }
//    std::cout << "\n";

    // Compute m_i_2 (we know it's perpendicular to e_i and m_i_1, so this shouldn't be too hard)
    float n[4][3];
    float n_equil[4][3];
    rod::cross_all(p, m, n);
    rod::cross_all(p_equil, m_equil, n_equil);
    
    // Compute unperturbed energy.
    // I could make this less verbose, but the explicit lookup table is a bit clearer about what's going on.
    // The basic idea is: if we're close to the 'edge' of the rod, don't compute energies for non-existent nodes! Because they are out of bounds!
    float bend_energy = 0;
    float stretch_energy = 0;
    float twist_energy = 0;
    
//    return 0;
//    std::cout << "start cutoff = " << start_cutoff << ", end cutoff = " << end_cutoff << "\n";
        
    if (start_cutoff == 0 && end_cutoff == 0){
        twist_energy += rod::get_twist_energy(material[i][1], m[i], m[im1], m_equil[i], m_equil[im1], p[im1], p[i], p_equil[im1], p_equil[i]);
        twist_energy += rod::get_twist_energy(material[ip1][1], m[ip1], m[i], m_equil[ip1], m_equil[i], p[i], p[ip1], p_equil[i], p_equil[ip1]);
        twist_energy += rod::get_twist_energy(material[im1][1], m[im1], m[im2], m_equil[im1], m_equil[im2], p[im2], p[im1], p_equil[im2], p_equil[im1]);
    }
    else if (start_cutoff == 1 && end_cutoff == 0){
        twist_energy += rod::get_twist_energy(material[i][1], m[i], m[im1], m_equil[i], m_equil[im1], p[im1], p[i], p_equil[im1], p_equil[i]);
        twist_energy += rod::get_twist_energy(material[ip1][1], m[ip1], m[i], m_equil[ip1], m_equil[i], p[i], p[ip1], p_equil[i], p_equil[ip1]);
    }
    else if (start_cutoff == 2 && end_cutoff == 0){
        twist_energy += rod::get_twist_energy(material[ip1][1], m[ip1], m[i], m_equil[ip1], m_equil[i], p[i], p[ip1], p_equil[i], p_equil[ip1]);
    }
    else if (start_cutoff == 0 && end_cutoff == 1){
        twist_energy += rod::get_twist_energy(material[i][1], m[i], m[im1], m_equil[i], m_equil[im1], p[im1], p[i], p_equil[im1], p_equil[i]);
        twist_energy += rod::get_twist_energy(material[im1][1], m[im1], m[im2], m_equil[im1], m_equil[im2], p[im2], p[im1], p_equil[im2], p_equil[im1]);
    }    
    else if (start_cutoff == 0 && end_cutoff == 2){
        twist_energy += rod::get_twist_energy(material[im1][1], m[im1], m[im2], m_equil[im1], m_equil[im2], p[im2], p[im1], p_equil[im2], p_equil[im1]);
     }
    else{
        throw std::invalid_argument("Length of the rod must be larger than 3");
    }

    return twist_energy;
} 

int main(){
    rod::Rod test_rod("twisted_bent_rod.rod", 0);
    test_rod.load_header("twisted_bent_rod.rod");
    test_rod.load_contents("twisted_bent_rod.rod");
    test_rod.set_units();
    
    float* currm_data = test_rod.current_m;
    
//    rod::print_array("rod current_r", currm_data, 18);
    
    float perturbation_amount = (0.001*pow(10,-9))/mesoDimensions::length;
    
    float test_energy_up = twist_energy_perturbation_test(test_rod, y, 0, perturbation_amount);
    float test_energy_down = twist_energy_perturbation_test(test_rod, y, 0, perturbation_amount*-1);
//    float test_energy_neutral = twist_energy_perturbation_test(test_rod, y, 0, 0);
    //std::cout << "Perturbation amount: " << test_rod.perturbation_amount << "\n";
    std::cout << "Energy up: " << test_energy_up << "\n";
    std::cout << "Energy down: " << test_energy_down << "\n";
//    std::cout << "Energy neutral: " << test_energy_neutral << "\n";
    return 0;
}
