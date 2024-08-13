#include "rod_math_v9.h"

int main(){
    // testing rodriguez
    float e_i[3] = {1,0,0};
    float m_i[3] = {0,1,0};
    float rotated_mi[3];
    float theta = 1.570796327*2;
    rod::rodrigues_rotation(m_i, e_i, theta, rotated_mi);
    if (rotated_mi[1] < -0.99 && rotated_mi[1] > -1.01 && rotated_mi[2] < 0.001 && rotated_mi[2] > -0.001){
        return 0;
    }
    return 1;
} 
