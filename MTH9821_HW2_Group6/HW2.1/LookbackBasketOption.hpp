//
//  LookbackBasketOption.hpp
//  HW2.1
//
//  Created by 王明森 on 2022/09/10.
//

#ifndef LookbackBasketOption_hpp
#define LookbackBasketOption_hpp

#include <tuple>
#include <vector>

class LookbackBasketOption {
private:
    double S0_1_;   // spot price of the first asset
    double S0_2_;   // spot price of the second asset
    double K_;      // strike price
    double T_;      // maturity
    double sigma_1_;// volatility, first asset
    double sigma_2_;// volatility, first asset
    double r_;      // interest rate (constant & coumpounded continuously)
    double q_1_;    // dividend rate, first asset (constant & coumpounded continuously)
    double q_2_;    // dividend rate, second asset (constant & coumpounded continuously)
    double rho_;    // correlation of the price of the two assets
    std::size_t steps_; // Number of steps
    
    // Helper constants
    double deltat_;         // T / steps
    double w2_;             // sqrt(1 - rho^2), weight of the second standard normal variable
    double drift_1_;        // (r - q - sigma_1^2/2) * dt
    double drift_2_;        // (r - q - sigma_2^2/2) * dt
    double diffusion_1_;    // sigma_1 * sqrt(T)
    double diffusion_2_;    // sigma_2 * sqrt(T)
    double r_disc_;         // Discounter
    
public:
    LookbackBasketOption(double S0_1, double S0_2, double K, double T, double sigma_1, double sigma_2, double r, double q_1, double q_2, double rho, std::size_t steps);
    ~LookbackBasketOption() = default;
    
    double Price(const std::vector<double>& z1, const std::vector<double>& z2) const;  // Takes in a pair of independent standard normal variables
    
};


#endif /* LookbackBasketOption_hpp */
