// FDM: multiple classes for approximating SDES
// // taken out of the Advanced C++ course
// @ MTH9821 Homework1 Group6 

#include "Fdm.h"
#include <cmath>

/* FDM BASE*/
// Constructor
FdmBase::FdmBase(const std::shared_ptr<ISde>& sde, unsigned NT) : sde_(sde), NT_(NT), dt_(sde->Expiry() / double(NT)) {
    dt_sqrt_ = std::sqrt(dt_);
}

// Get / set underlying SDE
std::shared_ptr<ISde> FdmBase::SDE() const {
    return sde_;
}

void FdmBase::SDE(const std::shared_ptr<ISde>& sde) {
    sde_ = sde;
}

/* EULER */
// One-factor explicit Euler
Euler::Euler(const std::shared_ptr<ISde>& sde, unsigned NT) : FdmBase(sde, NT) {}

double Euler::advance(double x, unsigned step, double wiener_increment) const {
    double t = double(step) * dt_;  // Find current time
    
    // x(t+dt) = x(t) + mu(x, t) * dt + sigma(x, t) * dB
    return x + sde_->Convection(x, t) * dt_ + sde_->Diffusion(x, t) * dt_sqrt_ * wiener_increment;
}

/* Multiplicative */
// Multiplicative formula
Multiplicative::Multiplicative(const std::shared_ptr<ISde>& sde, unsigned NT, double sig, double mu) : FdmBase(sde, NT), sig_(sig), mu_(mu) {}

double Multiplicative::advance(double x, unsigned step, double wiener_increment) const {
    double alpha = .5 * sig_ * sig_;
    
    // x(t+dt) = x(t) * exp((mu - sigma^2 / 2) * dt + sigma * dB)
    return x * std::exp((mu_ - alpha) * dt_ + sig_ * dt_sqrt_ * wiener_increment);
}

/* Exact */
// Compute exact value at t + dt
Exact::Exact(const std::shared_ptr<ISde>& sde, unsigned NT, double sig, double mu) : FdmBase(sde, NT), S0_(sde->S0()), sig_(sig), mu_(mu) {}

double Exact::advance(double x, unsigned step, double wiener_increment) const {
    double t = double(step) * dt_;  // Find current time
    
    double alpha = .5 * sig_ * sig_;
    
    // Compute the exact value at (t + dt)
    return S0_ * std::exp((mu_ - alpha) * (t + dt_) + sig_ * std::sqrt(t + dt_) * wiener_increment);
}

/* Milstein */
// Milstein method
Milstein::Milstein(const std::shared_ptr<ISde>& sde, unsigned NT) : FdmBase(sde, NT) {}

double Milstein::advance(double x, unsigned step, double wiener_increment) const {
    double t = double(step) * dt_;  // Find current time
    
    // x(t+dt) = x(t) + mu * dt + sigma * dB + .5 * mu * mu_derivative (dB^2 - dt)
    // https://en.wikipedia.org/wiki/Milstein_method#Intuitive_derivation
    return x + sde_->Convection(x, t) * dt_ + sde_->Diffusion(x, t) * dt_sqrt_ * wiener_increment + .5 * dt_ * sde_->Convection(x, t) * sde_->DiffusionDerivative(x, t) * (wiener_increment * wiener_increment - 1.);
}

/* Predictor-Corrector */
// Predictor-corrector method
PredictorCorrector::PredictorCorrector(const std::shared_ptr<ISde>& sde, unsigned NT, double a, double b) : FdmBase(sde, NT), a_(a), b_(b) {}

double PredictorCorrector::advance(double x, unsigned step, double wiener_increment) const {
    double t = double(step) * dt_;  // Find current time
    
    // https://en.wikipedia.org/wiki/Predictor–corrector_method
    
    // Euler for predictor
    double VMid = x + sde_->Convection(x, t) * dt_ + sde_->Diffusion(x, t) * dt_sqrt_ * wiener_increment;
    
    // Modified double trapezoidal rule
    double convection_double_term = (a_ * sde_->Convection(VMid, t + dt_) + ((1. - a_) * sde_->Convection(x, t))) * dt_;
    double diffusion_double_term = (b_ * sde_->Diffusion(VMid, t + dt_) + ((1. - b_) * sde_->Diffusion(x, t))) * dt_sqrt_ * wiener_increment;
    
    return x + convection_double_term + diffusion_double_term;
}

/* Modified Predictor-Corrector */
// Modified predictor-Corrector method
ModifiedPredictorCorrector::ModifiedPredictorCorrector(const std::shared_ptr<ISde>& sde, unsigned NT, double a, double b) : FdmBase(sde, NT), a_(a), b_(b) {}

double ModifiedPredictorCorrector::advance(double x, unsigned step, double wiener_increment) const {
    double t = double(step) * dt_;  // Find current time
    
    // Predictor-corrector using adjusted convection and trapezoidal rule
    // https://en.wikipedia.org/wiki/Predictor–corrector_method
    
    // Euler for predictor
    double VMid = x + sde_->Convection(x, t) * dt_ + sde_->Diffusion(x, t) * dt_sqrt_ * wiener_increment;
    
    // Modified double trapezoidal rule
    double convection_double_term = (a_ * sde_->ConvectionCorrected(VMid, t + dt_, b_) + ((1. - a_) * sde_->ConvectionCorrected(x, t, b_))) * dt_;  // Use corrected convection instead of naive convection
    double diffusion_double_term = (b_ * sde_->Diffusion(VMid, t + dt_) + ((1. - b_) * sde_->Diffusion(x, t))) * dt_sqrt_ * wiener_increment;
    
    return x + convection_double_term + diffusion_double_term;
}

