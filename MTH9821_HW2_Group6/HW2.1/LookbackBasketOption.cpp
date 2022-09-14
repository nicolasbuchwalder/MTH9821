//
//  LookbackBasketOption.cpp
//  HW2.1
//
//  Created by 王明森 on 2022/09/10.
//

#include "LookbackBasketOption.hpp"
#include <cmath>
#include <algorithm>

LookbackBasketOption::LookbackBasketOption(double S0_1, double S0_2, double K, double T, double sigma_1, double sigma_2, double r, double q_1, double q_2, double rho, std::size_t steps) : S0_1_(S0_1), S0_2_(S0_2), K_(K), T_(T), sigma_1_(sigma_1), sigma_2_(sigma_2), r_(r), q_1_(q_1), q_2_(q_2), rho_(rho), steps_(steps),
    deltat_(T / steps), w2_(std::sqrt(1. - rho * rho)),
    drift_1_((r - q_1 - sigma_1 * sigma_1 / 2.) * deltat_),
    drift_2_((r - q_2 - sigma_2 * sigma_2 / 2.) * deltat_),
    diffusion_1_(sigma_1 * std::sqrt(deltat_)),
    diffusion_2_(sigma_2 * std::sqrt(deltat_)),
    r_disc_(std::exp(-r * T)) {}

double LookbackBasketOption::Price(const std::vector<double>& z1, const std::vector<double>& z2) const {
    // The original spot price
    double S1 = S0_1_;
    double S2 = S0_2_;
    
    // S1 + S2
    double running_max = S1 + S2;
    
    auto z1it = z1.cbegin();
    auto z2it = z2.cbegin();
    
    for (std::size_t j = 0; j < steps_; j++) {
        S1 *= std::exp(drift_1_ + diffusion_1_ * (*z1it));
        S2 *= std::exp(drift_2_ + diffusion_2_ * ((*z1it) * rho_ + (*z2it) * w2_));
        
        if (running_max < (S1 + S2)) {
            running_max = (S1 + S2);
        }
        
        z1it++;
        z2it++;
    }
    
    return r_disc_ * std::max(running_max - K_, 0.);
}
