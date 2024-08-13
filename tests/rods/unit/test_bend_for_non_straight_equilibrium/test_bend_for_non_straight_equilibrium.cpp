#include "rod_structure.h"
#include <iostream>
#include <string>

int main(){
    float mathematica_result = 0.08;
    
    float e_im1[3] = {50, 0, 0};
    float e_i[3] = {50, 0, 0};
    float e_im1_bar[3] = {50, 0, 0};
    float e_i_bar[3] = {0, 50, 0};

    float m_im1[3] = {0, 0, 1};
    float m_im1_bar[3] = {0, 0, 1};
    float m_i[3] = {0, 0, 1};
    float m_i_bar[3] = {0, 0, 1};
    
    float m_im1_2[3];
    float m_im1_bar_2[3];
    float m_i_2[3];
    float m_i_bar_2[3];
    
    rod::cross_product(m_im1, e_im1, m_im1_2);
    rod::cross_product(m_im1_bar, e_im1_bar, m_im1_bar_2);
    rod::cross_product(m_i, e_i, m_i_2);
    rod::cross_product(m_i_bar, e_i_bar, m_i_bar_2);
    
    rod::normalize(m_im1_2, m_im1_2);
    rod::normalize(m_im1_bar_2, m_im1_bar_2);
    rod::normalize(m_i_2, m_i_2);
    rod::normalize(m_i_bar_2, m_i_bar_2);
    
    float B[4] = {1,0,0,1};

    //float kbi[3];
    //float kbibar[3];
    
    float computed_energy = rod::get_bend_energy_from_p(
        e_im1,
        e_i,
        e_im1_bar,
        e_i_bar,
        m_im1_2,
        m_im1,
        m_im1_bar_2,
        m_im1_bar,
        m_i_2,
        m_i,
        m_i_bar_2,
        m_i_bar,
        B,
        B);
        
    if (computed_energy > mathematica_result - 0.00001 && computed_energy < mathematica_result + 0.00001){
        return 0;
    }
    else{
        std::cout << "Computed energy: " << computed_energy << "\n";
        std::cout << "Analytical energy: " << mathematica_result << "\n";
        return 1;
    }
    
    //get_bend_energy_from_e(
    //get_kb_i(eim1, ei, kbi);
    //get_kb_i(eim1bar, eibar, kbibar);
    
}
