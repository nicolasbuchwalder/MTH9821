//
//  main.cpp
//  HW2.1
//
//  Created by 王明森 on 2022/09/09.
//

#include <iostream>
#include "OneStepMC.hpp"
#include "MultiStepMC.hpp"

int main(int argc, const char * argv[]) {
    
    // Initialize the simulator (with pre-generated uniform RVs)
    OneStepMC simulator((1 << 9) * 10000);
    
    // set vanilla option data
    EuropeanOption option(56., 54., .75, .27, .02, 0.);

    // simulate with control variate
    simulator.ControlVariate(option);

    // simulate with antithetic variables
    simulator.AntitheticVariable(option);

    // simulate with moment matching
    simulator.MomentMatching(option);

    // simulate with moment matching and control variates
    simulator.MMCV(option);
    
    // set basket option data
    BasketOption basket(26., 29., 50., .5, .31, .21, .025, 0., 0., .3);

    // find price using mc
    simulator.Price(basket);
    
    
    // Initialize the path-dependent simulator
    MultiStepMC path_simulator(150, 50 * (1 << 9), 2);

    // Set lookback basket option data
    LookbackBasketOption lookback_basket(26., 29., 50., .5, .31, .21, .025, 0., 0., .3, 150);

    // Find price using MC
    path_simulator.Price(lookback_basket);

    
    return 0;
}
