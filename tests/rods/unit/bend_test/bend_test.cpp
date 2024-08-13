#include "rod_math_v9.h"

float get_bend_for_test(float x2_current[3]){
    float x1[3] = {0,0,0};
    float x2[3] = {0,1,0};
    float x3[3] = {0,2,0};
    float eim1bar[3];
    float eibar[3];
    rod::get_p_i(x1, x2, eim1bar);
    rod::get_p_i(x2, x3, eibar);
    float mibar[3];
    float mim1bar[3];
    rod::cross_product(eibar, eim1bar, mibar);
    rod::cross_product(eibar, eim1bar, mim1bar);
    float mi2bar[3];
    float mim2bar[3];
    rod::cross_product(eibar, mibar, mi2bar);
    rod::cross_product(eim1bar, mim1bar, mim2bar);
    
    float x1_current[3] = {0,0,0};
//    float x2_current[3] = {0,1,-1};
    float x3_current[3] = {0,2,0};
    float eim1[3];
    float ei[3];
    rod::get_p_i(x1_current, x2_current, eim1);
    rod::get_p_i(x2_current, x3_current, ei);
    float mi[3];
    float mim1[3];
    rod::cross_product(ei, eim1, mi);
    rod::cross_product(ei, eim1, mim1);
    float mi2[3];
    float mim2[3];
    rod::cross_product(ei, mi, mi2);
    rod::cross_product(eim1, mim1, mim2);
       
    float x2_perturbed_x [3] = {0.05,1,-1};
    float x2_perturbed_y [3] = {0,1.05,-1};
    float x2_perturbed_z [3] = {0,1,-0.95};

    float e_i_x2_x[3];
    float e_i_x2_y[3];
    float e_i_x2_z[3];
    
    float e_im1_x2_x[3];
    float e_im1_x2_y[3];
    float e_im1_x2_z[3];

    rod::get_p_i(x1, x2_perturbed_x, e_im1_x2_x);
    rod::get_p_i(x1, x2_perturbed_y, e_im1_x2_y);
    rod::get_p_i(x1, x2_perturbed_z, e_im1_x2_z);
    
    rod::get_p_i(x2_perturbed_x, x3, e_i_x2_x);
    rod::get_p_i(x2_perturbed_y, x3, e_i_x2_y);
    rod::get_p_i(x2_perturbed_z, x3, e_i_x2_z);
    
    rod::get_p_i(x2, x3, eibar);
    
    float B_i_equil[4] = {1,0,0,1};
    
    float bend_energy_x = rod::get_bend_energy_from_p(e_im1_x2_x, e_i_x2_x, eim1bar, eibar, mim2, mim1, mim2bar, mim1, mi2, mi, mi2bar, mibar, B_i_equil, B_i_equil);
    float bend_energy_y = rod::get_bend_energy_from_p(e_im1_x2_y, e_i_x2_y, eim1bar, eibar, mim2, mim1, mim2bar, mim1, mi2, mi, mi2bar, mibar, B_i_equil, B_i_equil);
    float bend_energy_z = rod::get_bend_energy_from_p(e_im1_x2_z, e_i_x2_z, eim1bar, eibar, mim2, mim1, mim2bar, mim1, mi2, mi, mi2bar, mibar, B_i_equil, B_i_equil);

    return bend_energy_x+bend_energy_y+bend_energy_z;
}

int main(){
    float x2_current_90[3] = {0,1,-1};
    float x2_current_45[3] = {0,1,-0.5};
    float energy_90 = get_bend_for_test(x2_current_90);
    float energy_45 = get_bend_for_test(x2_current_45);
    if (energy_90 > energy_45 && energy_90 < 50 && energy_45 < 50){
        return 0;
    }
    return 1;
}
