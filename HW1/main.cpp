//
//  main.cpp
//  HW1
//
//  Created by 王明森 on 2022/08/30.
//

#include <iostream>
#include <vector>
#include "OneStepMC.hpp"
#include "PseudoRNG.hpp"
#include <numeric>

int main(int argc, const char * argv[]) {
    // S0, K, T, sig, r, q
    EuroOption option(41., 42., .75, .25, .03, .01);

    OneStepMC mc;

//    mc.DirectSimulation(option);
    
    double deltaS = 0.01;
    EuroOption lower(41. - deltaS, 42., .75, .25, .03, .01);
    EuroOption upper(41. + deltaS, 42., .75, .25, .03, .01);

    mc.FDMSimulation(option, lower, upper, deltaS);

    return 0;
}
