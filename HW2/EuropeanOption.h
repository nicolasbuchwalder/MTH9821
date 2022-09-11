// EuropeanOption: class that represents an European option (call or put)
// @ MTH9821 Homework2 Group6 

#ifndef EUROPEANOPTION_H
#define EUROPEANOPTION_H

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

    // useful values
    double d1_;
    double d2_;
    double Zd1_;
    double Zd2_;
    double Nd1_;
    double Nd2_;
    double r_disc_; // discount by interest rate: exp(-r(T-t))
    double q_disc_; // discount by dividend rate: exp(-q(T-t))

    // helper functions
    double Z(double t) const;   // PDF of a standard normal variable
    double Phi(double t) const; // CDF of a standard normal variable

public:
    // value constructor with characteristics of option
    EuropeanOption(double S0, double K, double T, double sigma, double r, double q);

    ~EuropeanOption() = default;

    // function that generates all values of prices and greeks at maturity for a given z 
    std::array<double, 6> Outcome(double z) const;


    double CallFromNormal(double z) const;
    double PutFromNormal(double z) const;

    double CallFromStock(double S) const;
    double PutFromStock(double S) const;

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
};

#endif /* EUROPEANOPTION_H */
