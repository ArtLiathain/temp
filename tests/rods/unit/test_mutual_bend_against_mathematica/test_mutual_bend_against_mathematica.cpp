#include "rod_math_v9.h"

int main(){
    //float mathematica_result = 059.6893; //from previous iteration of the algorithm
    float mathematica_result = 45.99122923;
    float bbar[4] = {1,0,0,1};
    float eim1[3] = {0,0, 0.008};
    float ei[3] = {0, 0.00841471, 0.00540302};
    float eim1bar[3] = {0,0,0.01};
    float eibar[3] = {0,0,0.01};
    float mim12[3] = {1,0,0};
    float mim1[3] = {0, 1, 0};
    float mim12bar[3] = {1,0,0};
    float mim1bar[3] = {0,1,0};
    float mi2[3] = {1,0,0};
    float mi[3] = {0, -0.540302, -0.841471};
    float mi2bar[3] = {1,0,0};
    float mibar[3] = {0,1,0};
    //float kbi[3];
    //float kbibar[3];
    
    float computed_energy = rod::get_bend_energy_mutual_parallel_transport(
        eim1,
        ei,
        eim1bar,
        eibar,
        mim12,
        mim1,
        mim12bar,
        mim1bar,
        mi2,
        mi,
        mi2bar,
        mibar,
        bbar,
        bbar);
        
    if (computed_energy > mathematica_result - 0.01 && computed_energy < mathematica_result + 0.01){
        std::cout << "Mutual bend against mathematica: passed\n";
        return 0;
    } else {
        std::cout << "Computed energy: " << computed_energy << "\n";
        std::cout << "Analytical energy: " << mathematica_result << "\n";
        return 1;
    }
    
    //get_bend_energy_from_e(
    //get_kb_i(eim1, ei, kbi);
    //get_kb_i(eim1bar, eibar, kbibar);
    
}
