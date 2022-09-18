// EuropeanOption: class that represents an European option (call or put)
// @ MTH9821 Homework2 Group6
#ifndef EUROPEANOPTION_hpp
#define EUROPEANOPTION_hpp

#include <array>
#include <vector>

class EuropeanOption {
private:
    double S0_;     // spot price
    double K_;      // strike price
    double T_;      // maturity
    double sigma_;  // volatility (constant)
    double r_;      // interest rate (constant & coumpounded continuously)
    double q_;      // dividend rate (constant & coumpounded continuously)

    double t_;

    // Values for the Black-Scholes Model
    double d1_;
    double d2_;
    double Zd1_;
    double Zd2_;
    double Nd1_;
    double Nd2_;
    double r_disc_; // Discount by interest rate: exp(-r(T-t))
    double q_disc_; // Discount by dividend rate: exp(-q(T-t))
    
    // helper functions
    double Z(double t) const;   // PDF of a standard normal variable
    double Phi(double t) const; // CDF of a standard normal variable

public:
    // value constructor with characteristics of option
    EuropeanOption(double S0, double K, double T, double sigma, double r, double q);

    ~EuropeanOption() = default;

    // function that generates all values of prices and greeks at maturity for a given z
//    std::array<double, 6> OptionGreeks(double z) const;

    // Spot price and the price of the put
    std::array<double, 2> SpotAndPut(double z) const;

    double Call(double z) const;
    double Put(double z) const;

    // BLACK-SCHOLES
    
    // price of both call and put options with given specificities
    double Call() const;
    double Put() const;

    // Greeks
    double DeltaCall() const;
    double DeltaPut() const;
    double GammaCall() const;
    double GammaPut() const;
    double VegaCall() const;
    double VegaPut() const;
    
    friend class OneStepMC;
};

#endif /* EUROPEANOPTION_hpp */
