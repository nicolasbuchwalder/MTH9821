// ISde: class that evaluates Stochastic differential equations
// taken out of the Advanced C++ course
// @ MTH9821 Homework1 Group6 

#include "Sde.h"
#include <cmath>

/* Generic SDE */
ISde::ISde(double S0, double expiry) : S0_(S0), expiry_(expiry) {}

// Get / set expiry time
double ISde::Expiry() const { return expiry_; }
void ISde::Expiry(double expiry) { expiry_ = expiry; }

// Get / set initial condition
double ISde::S0() const { return S0_; }
void ISde::S0(double s0) { S0_ = s0; }

/* GEOMETRIC BROWNIAN MOTION */
// Constructor
GBM::GBM(double riskless_rate, double volatility, double dividend_rate, double initial_condition, double expiry) :
    ISde(initial_condition, expiry), r_(riskless_rate), vol_(volatility), d_(dividend_rate) {}

// Coefficients
double GBM::Convection(double x, double t) const {
    return (r_ - d_) * x;
}
double GBM::Diffusion(double x, double t) const {
    return vol_ * x;
}
double GBM::ConvectionCorrected(double x, double t, double B) const {
    return (this->Convection(x, t) - B * Diffusion(x, t) * DiffusionDerivative(x, t));
}
double GBM::DiffusionDerivative(double x, double t) const {
    return vol_;
}

/* CONSTANT ELASTICITY OF VARIANCE */
// Constructor
CEV::CEV(double riskless_rate, double volatility, double dividend_rate, double initial_condition, double expiry, double beta) :
ISde(initial_condition, expiry), r_(riskless_rate), vol_(volatility), d_(dividend_rate), beta_(beta) {}

// Coefficients
double CEV::Convection(double x, double t) const {
    return (r_ - d_) * x;
}
double CEV::Diffusion(double x, double t) const {
    return vol_ * std::pow(x, beta_);
}
double CEV::ConvectionCorrected(double x, double t, double B) const {
    return (this->Convection(x, t) - B * Diffusion(x, t) * DiffusionDerivative(x, t));
}
double CEV::DiffusionDerivative(double x, double t) const {
    return vol_ * beta_ * std::pow(x, beta_ - 1.);
}

