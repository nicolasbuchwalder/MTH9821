//
//  OneStepMC.hpp
//  HW1
//
//  Created by 王明森 on 2022/08/30.
//

#ifndef OneStepMC_hpp
#define OneStepMC_hpp

#include <array>
#include <vector>

class EuroOption {
private:
    double S0_;     // Spot price
    double K_;      // Strike price
    double T_;      // Maturity
    double sigma_;  // Volatility
    double r_;      // Const interest rate
    double q_;      // Dividend rate
    
    double t_;
    
    // Intermediate values
    double d1_;
    double d2_;
    double Zd1_;
    double Zd2_;
    double Nd1_;
    double Nd2_;
    double q_disc_; // Discount by dividend rate: exp(-q(T-t))
    double r_disc_; // Discount by interest rate: exp(-r(T-t))
    
    double Z(double t) const;   // PDF of a standard normal variable
    double Phi(double t) const; // CDF of a standard normal variable
    
public:
    EuroOption(double S0, double K, double T, double sigma, double r, double q);
    ~EuroOption() = default;
    
    // Generate from standard normal distribution
    std::array<double, 6> Outcome(double z) const;
    double Call(double z) const;
    
    // BLACK-SCHOLES
    // Price
    double Call() const;
    double Put() const;
    
    // Greeks
    double DeltaCall() const;
    double DeltaPut() const;
    double GammaCall() const;
    double GammaPut() const;
    double VegaCall() const;
    double VegaPut() const;
};

class OneStepMC {
private:
    std::vector<double> z_;
    
public:
//    OneStepMC(const EuroOption& option);
    OneStepMC() = default;
    ~OneStepMC() = default;
    
    void DirectSimulation(const EuroOption& option);
    void FDMSimulation(const EuroOption& option, const EuroOption& lower, const EuroOption& upper, double deltaS);
};

#endif /* OneStepMC_hpp */
