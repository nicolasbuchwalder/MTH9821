//
//  BasketOption.cpp
//  HW2.1
//
//  Created by 王明森 on 2022/09/10.
//

#include "BasketOption.hpp"
#include <cmath>
#include <algorithm>

BasketOption::BasketOption(double S0_1, double S0_2, double K, double T, double sigma_1, double sigma_2, double r, double q_1, double q_2, double rho) : S0_1_(S0_1), S0_2_(S0_2), K_(K), T_(T), sigma_1_(sigma_1), sigma_2_(sigma_2), r_(r), q_1_(q_1), q_2_(q_2), rho_(rho),
    w2_(std::sqrt(1. - rho * rho)),
    drift_1_((r - q_1 - sigma_1 * sigma_1 / 2.) * T),
    drift_2_((r - q_2 - sigma_2 * sigma_2 / 2.) * T),
    diffusion_1_(sigma_1 * std::sqrt(T)),
    diffusion_2_(sigma_2 * std::sqrt(T)),
    r_disc_(std::exp(-r * T)) {}

double BasketOption::Price(double z1, double z2) const {
    double S1 = S0_1_ * std::exp(drift_1_ + diffusion_1_ * z1);
    double S2 = S0_2_ * std::exp(drift_2_ + diffusion_2_ * (rho_ * z1 + w2_ * z2));
    return r_disc_ * std::max(S1 + S2 - K_, 0.);
}
