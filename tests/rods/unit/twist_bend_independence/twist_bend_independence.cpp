#include "rod_math_v9.h"

int main(){
    float m_i[3] = {0,1,0};
    float m_im1[3] = {0,1,0};
    float m_i_equil[3] = {0,1,0};
    float m_im1_equil[3] = {0,1,0};
    float e_im1[3] = {1,0,0};
    float e_i[3] = {1,0,0};
    float e_im1_equil[3] = {1,0,0};
    float e_i_equil[3] = {1,0,0};
    rod::rodrigues_rotation(m_i, e_i, 4.2, m_i);
    float twist = rod::get_twist_energy(1, m_i, m_im1, m_i_equil, m_im1_equil, e_im1, e_i, e_im1_equil, e_i_equil);
    
    float m_i_2[3];
    rod::cross_product(m_i, e_i, m_i_2);

    float m_i_2_equil[3];
    rod::cross_product(m_i_equil, e_i_equil, m_i_2_equil);
    
    float m_im1_2[3];
    rod::cross_product(m_im1, e_im1, m_im1_2);

    float m_im1_2_equil[3];
    rod::cross_product(m_im1_equil, e_im1_equil, m_im1_2_equil);
    
    float B_equil[4] = {1,0,0,1};
    
    float bend = rod::get_bend_energy_from_p(e_im1, e_i, e_im1_equil, e_i_equil, m_im1_2, m_im1, m_im1_2_equil,
        m_im1_equil,
        m_i_2,
        m_i,
        m_i_2_equil,
        m_i_equil,
        B_equil,
        B_equil);

    if (twist > 0 && bend > -0.001 && bend < 0.001){
        return 0;
    }

//    std::cout << "m_i rotated = " << m_i[0] << ", " << m_i[1] << ", " << m_i[2] << "\n";
//    std::cout << "twist emergy (should be >0)" << twist << "\n";
//    std::cout << "bend energy (should be ~0)" << bend << "\n";

    return 1;
} 
