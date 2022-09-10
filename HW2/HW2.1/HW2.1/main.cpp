//
//  main.cpp
//  HW2.1
//
//  Created by 王明森 on 2022/09/09.
//

#include <iostream>
#include "EuropeanOption.hpp"
#include "OneStepMC.hpp"

int main(int argc, const char * argv[]) {
    
    // Initialize the simulator (with pre-generated uniform RVs)
    OneStepMC simulator((1 << 9) * 10000);
    
    // Set option data
    EuropeanOption option(56., 54., .75, .27, .02, 0.);
    
    // Simulate with Control Variate
    simulator.ControlVariate(option);
    
    // Simulate with Antithetic Variables
    simulator.AntitheticVariable(option);
    
    // Simulate with Moment Matching
    simulator.MomentMatching(option);
    
    // Simulate with Moment Matching and Control Variates
    simulator.MMCV(option);
    
    return 0;
}
