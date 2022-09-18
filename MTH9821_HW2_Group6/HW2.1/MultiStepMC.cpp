//
//  MultiStepMC.cpp
//  HW2.1
//
//  Created by 王明森 on 2022/09/10.
//

#include "MultiStepMC.hpp"
#include "LCG.hpp"
#include "BoxMuller.hpp"
#include <iomanip>
#include <iostream>

MultiStepMC::MultiStepMC(std::size_t steps, std::size_t paths, std::size_t multiplier) {
    // Generate samples using Box Muller
    LCG::reseed(1); // Reset LCE
    z_.resize(paths * multiplier);
    
    for (auto zgen_it = z_.begin(); zgen_it != z_.end(); zgen_it++) {
        zgen_it->resize(steps);
        
        auto pathit = zgen_it->begin();
        
        // Just in case there are an odd number of steps
        if (steps % 2) {
            double z1, z2;
            std::tie(z1, z2) = BoxMuller::standard_normal_pair();
            *(pathit++) = z1;
        }
        
        while (pathit != zgen_it->end()) {
            double z1, z2;
            std::tie(z1, z2) = BoxMuller::standard_normal_pair();
            *(pathit++) = z1;
            *(pathit++) = z2;
        }
    }
}

void MultiStepMC::Price(const LookbackBasketOption& option) const {
    // Prepare for Monte Carlo
    double V_sum = 0.;
    
    // BEGINNING AND END OF SUBSEQUENCES
    std::vector<std::size_t> begins({ 0, 1, 2, 4, 8, 16, 32, 64, 128, 256 });
    std::vector<std::size_t> ends({ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512 });
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "======= LOOKBACK BASKET =======" << std::endl;
    
    auto zit = z_.cbegin();
    
    for (std::size_t loop = 0; loop < begins.size(); loop++) {

        // GENERATE PRICES AND GREEKS FOR EVERY ITERATION
        for (std::size_t i = begins[loop] * 50; i < ends[loop] * 50; i++) {
            
            double V = option.Price(*(zit++), *(zit++));
            
            V_sum += V;
            
        }

        double N = ends[loop] * 50;
        
        // FIND THE CURRENT MEAN (MONTE CARLO ESTIMATE)
        double V_hat = V_sum / N;

        // OUTPUT RESULTS

        std::cout << std::setw(5) << ends[loop] * 50 <<": ";

        // call results
        std::cout << std::setw(8) << V_hat << std::endl;
    }
    
    std::cout << std::endl;
}

