//
//  Mediator.cpp
//  Final_Project
//
//  Created by Wanderers' Library on 2022/2/26.
//

#include "Mediator.h"
#include <functional>

Mediator::Mediator(const std::shared_ptr<ISde>& sde, const std::shared_ptr<FdmBase>& fdm, unsigned long NSim, unsigned NT) : sde_(sde), fdm_(fdm), NSim_(NSim), NT_(NT) {}

void Mediator::AddMis(const MIS& mis) {
    progression_signal_.connect(mis);
}

void Mediator::AddPricer(const std::shared_ptr<Pricer>& pricer) {
    // Connect signals for pricing, post-processing, and printing results.
    
    pricer_signal_.connect(std::bind(&Pricer::ProcessPath, pricer, std::placeholders::_1));
    
    stop_signal_.connect(std::bind(&Pricer::PostProcess, pricer));
    
    print_signal_.connect(std::bind(&Pricer::PrintResults, pricer));
}

void Mediator::simulate() {
    
    
    for (unsigned long sim = 1; sim <= NSim_; sim++) {
        // Draw a path for each simulation
        
        std::vector<double> path;   // Simulated path
        
        path.reserve(NT_ + 1);         // The path's length will be time steps + 1 (the initial condition)
        path.push_back(sde_->S0());    // The path's beginning is determined by the initial condition
        
        if ((sim / 10000 * 10000) == sim) progression_signal_(sim);
        
        for (unsigned time = 0; time < NT_; time++) {
            // Find the new value for every time step
            double VNew = fdm_->advance(path.back(), time, ITM::standard_normal());    // Advance time from t to t+1
            path.push_back(VNew);
        }
        
        pricer_signal_(path);  // Send copies of the path to each pricer
    }
    
    stop_signal_(); // Tell pricers to stop and return results.
}

void Mediator::PrintResults() const {
    print_signal_();
}
