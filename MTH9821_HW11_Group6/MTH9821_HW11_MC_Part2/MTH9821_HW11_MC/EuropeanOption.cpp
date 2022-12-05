//
//  EuropeanOption.cpp
//  MonteCarlo1
//
//  Created by Aubree Li on 11/22/22.
//

#include "EuropeanOption.hpp"

#include <cmath>
#include <numbers>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

// value constructor with characteristics of option
EuropeanOption::EuropeanOption(double S0, double K, double T, double sigma, double r, double q)
    : S0_(S0), K_(K), T_(T), sigma_(sigma), r_(r), q_(q), t_(0.)
{
    double t = 0.;

    // calculating useful values for price and greeks calculations
    d1_ = (log(S0 / K) + (r - q + sigma * sigma / 2.) * (T - t)) / (sigma * std::sqrt(T - t));
    d2_ = d1_ - sigma * std::sqrt(T - t);

    Zd1_ = this->Z(d1_);
    Zd2_ = this->Z(d2_);
    Nd1_ = this->Phi(d1_);
    Nd2_ = this->Phi(d2_);

    q_disc_ = std::exp(-q * (T - t));
    r_disc_ = std::exp(-r * (T - t));
}

// function that generates all values of prices and greeks at maturity for a given z
std::array<double, 6> EuropeanOption::Outcome(double z) const {

    double sqrtT = std::sqrt(T_);                                                           // square root of time

    double S = S0_ * std::exp((r_ - q_ - sigma_ * sigma_ / 2.) * T_ + sigma_ * sqrtT * z);  // value of underlying

    double disc = std::exp(-r_ * T_);                                                       // discount factor

    double C = disc * std::max(S - K_, 0.);                                                 // value of call

    double DeltaC = ((S > K_) ? (disc * S / S0_) : 0.);                                     // delta of call

    double VegaC = ((S > K_) ? (S * disc * (-sigma_ * T_ + sqrtT * z)) : 0.);               // vega of call

    double P = disc * std::max(K_ - S, 0.);                                                 // value of put

    double DeltaP = ((S < K_) ? (-disc * S / S0_) : 0.);                                    // delta of put

    double VegaP = ((S < K_) ? (-S * disc * (-sigma_ * T_ + sqrtT * z)) : 0.);              // vega of put

    return std::array<double, 6>({ C, DeltaC, VegaC, P, DeltaP, VegaP });
}


double EuropeanOption::Call(double z) const {
    double sqrtT = std::sqrt(T_);                                                           // square root of time

    double S = S0_ * std::exp((r_ - q_ - sigma_ * sigma_ / 2.) * T_ + sigma_ * sqrtT * z);  // value of underlying

    double disc = std::exp(-r_ * T_);                                                       // discount factor

    double C = disc * std::max(S - K_, 0.);                                                 // value of call

    return C;
}

double EuropeanOption::Put(double z) const {
    double sqrtT = std::sqrt(T_);                                                           // square root of time

    double S = S0_ * std::exp((r_ - q_ - sigma_ * sigma_ / 2.) * T_ + sigma_ * sqrtT * z);  // value of underlying

    double disc = std::exp(-r_ * T_);                                                       // discount factor

    double P = disc * std::max(K_ - S, 0.);                                                // value of call

    return P;
}

// helper functions

// PDF of a z for random normal
double EuropeanOption::Z(double t) const {
    return std::exp(-t * t / 2.) / std::sqrt(2. * M_PI);
}

// CDF of a z for random normal
double EuropeanOption::Phi(double t) const {
    return std::erfc(-t / std::sqrt(2.)) / 2.;
}



// analytical solutions for prices and greeks

double EuropeanOption::Call() const {
    return S0_ * q_disc_ * Nd1_ - K_ * r_disc_ * Nd2_;
}

double EuropeanOption::Put() const {
    return -S0_ * q_disc_ * (1. - Nd1_) + K_ * r_disc_ * (1. - Nd2_);
}

double EuropeanOption::DeltaCall() const {
    return q_disc_ * Nd1_;
}
double EuropeanOption::DeltaPut() const {
    return -q_disc_ * (1. - Nd1_);
}

double EuropeanOption::GammaCall() const {
    return q_disc_ / (S0_ * sigma_ * std::sqrt(T_ - t_)) * Zd1_;
}
double EuropeanOption::GammaPut() const {
    return this->GammaCall();
}

double EuropeanOption::VegaCall() const {
    return S0_ * q_disc_ * std::sqrt(T_ - t_) * Zd1_;
}
double EuropeanOption::VegaPut() const {
    return this->VegaCall();
}
