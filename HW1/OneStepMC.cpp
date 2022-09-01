//
//  OneStepMC.cpp
//  HW1
//
//  Created by 王明森 on 2022/08/30.
//

#include "OneStepMC.hpp"
#include <cmath>
#include "PseudoRNG.hpp"
#include <iostream>
#include <numbers>
#include <iomanip>


EuroOption::EuroOption(double S0, double K, double T, double sigma, double r, double q) : S0_(S0), K_(K), T_(T), sigma_(sigma), r_(r), q_(q), t_(0.) {
    double t = 0;
    
    // The intermediate values are for BS price / Greek calculation
    d1_ = (log(S0 / K) + (r - q + sigma * sigma / 2.) * (T - t)) / (sigma * std::sqrt(T-t));
    d2_ = d1_ - sigma * std::sqrt(T-t);
    
    Zd1_ = this->Z(d1_);
    Zd2_ = this->Z(d2_);
    Nd1_ = this->Phi(d1_);
    Nd2_ = this->Phi(d2_);
    
    q_disc_ = std::exp(-q * (T-t));
    r_disc_ = std::exp(-r * (T-t));
}
std::array<double, 6> EuroOption::Outcome(double z) const {
    double sqrtT = std::sqrt(T_);
    
    double S = S0_ * std::exp((r_ - q_ - sigma_ * sigma_ / 2.) * T_ + sigma_ * sqrtT * z);
    
    double disc = std::exp(-r_ * T_);
    
    double C = disc * std::max(S - K_, 0.);
    
    double DeltaC = ((S > K_) ? (disc * S / S0_) : 0.);
    
    double VegaC = ((S > K_) ? (S * disc * (-sigma_ * T_ + sqrtT * z)) : 0.);
    
    double P = disc * std::max(K_ - S, 0.);
    
    double DeltaP = ((S < K_) ? (-disc * S / S0_) : 0.);
    
    double VegaP = ((S < K_) ? (-S * disc * (-sigma_ * T_ + sqrtT * z)) : 0.);
    
    return std::array<double, 6>({C, DeltaC, VegaC, P, DeltaP, VegaP});
}

double EuroOption::Call(double z) const {
    double sqrtT = std::sqrt(T_);
    
    double S = S0_ * std::exp((r_ - q_ - sigma_ * sigma_ / 2.) * T_ + sigma_ * sqrtT * z);
    
    double disc = std::exp(-r_ * T_);
    
    double C = disc * std::max(S - K_, 0.);
    
    return C;
}

double EuroOption::Phi(double t) const {
    return std::erfc(-t / std::sqrt(2.)) / 2.;
}

double EuroOption::Z(double t) const {
    return std::exp(-t * t / 2.) / std::sqrt(2. * std::numbers::pi);
}



double EuroOption::Call() const {
    return S0_ * q_disc_ * Nd1_ - K_ * r_disc_ * Nd2_;
}

double EuroOption::Put() const {
    return -S0_ * q_disc_ * (1. - Nd1_) + K_ * r_disc_ * (1. - Nd2_);
}

double EuroOption::DeltaCall() const {
    return q_disc_ * Nd1_;
}
double EuroOption::DeltaPut() const {
    return -q_disc_ * (1. - Nd1_);
}

double EuroOption::GammaCall() const {
    return q_disc_ / (S0_ * sigma_ * std::sqrt(T_ - t_)) * Zd1_;
}
double EuroOption::GammaPut() const {
    return this->GammaCall();
}



double EuroOption::VegaCall() const {
    return S0_ * q_disc_ * std::sqrt(T_ - t_) * Zd1_;
}

double EuroOption::VegaPut() const {
    return this->VegaCall();
}


//OneStepMC::OneStepMC(const EuroOption& option) : option_(option) {}


void OneStepMC::DirectSimulation(const EuroOption& option) {
    
    // GENERATE SAMPLES
    z_.resize(512 * 10000);
    
    std::generate(z_.begin(), z_.end(), BSM::standard_normal);
    
    auto zit = z_.cbegin();
    
    // PREPARE FOR MONTE CARLO
    std::array<double, 6> sum_of_outcome;
    std::fill(sum_of_outcome.begin(), sum_of_outcome.end(), 0.);
    
    // BEGINNING AND END OF SUBSEQUENCES
    std::vector<std::size_t> begins({0, 1, 2, 4, 8, 16, 32, 64, 128, 256});
    std::vector<std::size_t> ends({1, 2, 4, 8, 16, 32, 64, 128, 256, 512});
    
    std::cout << std::fixed << std::setprecision(6);
    
    for (std::size_t loop = 0; loop < begins.size(); loop++) {
        
        // GENERATE PRICES AND GREEKS FOR EVERY ITERATION
        for (std::size_t i = begins[loop] * 10000; i < ends[loop] * 10000; i++) {
            std::array<double, 6> out = option.Outcome(*(zit++));
            
            for (std::size_t j = 0; j < 6; j++) {
                sum_of_outcome[j] += out[j];
            }
            
        }
    
        // FIND THE CURRENT MEAN (MONTE CARLO ESTIMATE)
        std::array<double, 6> mean;
        std::transform(sum_of_outcome.cbegin(), sum_of_outcome.cend(), mean.begin(), [&](double sum)->double {
            return sum / (ends[loop] * 10000.);
        });
        
        // OUTPUT RESULTS
        std::cout << "======= " << ends[loop] << " =======" << std::endl;
        std::cout << "C, error(C), Delta(C), error(Delta(C)), Vega(C), error(Vega(C))" << std::endl;
        std::cout << mean[0] << std::endl;
        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[0] - option.Call()) << std::endl;
        std::cout << mean[1] << std::endl;
        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[1] - option.DeltaCall()) << std::endl;
        std::cout << mean[2] << std::endl;
        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[2] - option.VegaCall()) << std::endl;
        
        std::cout << "\nP, error(P), Delta(P), error(Delta(P)), Vega(P), error(Vega(P))" << std::endl;
        std::cout << mean[3] << std::endl;
        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[3] - option.Put()) << std::endl;
        std::cout << mean[4] << std::endl;
        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[4] - option.DeltaPut()) << std::endl;
        std::cout << mean[5] << std::endl;
        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[5] - option.VegaPut()) << std::endl;
    }
    
}

void OneStepMC::FDMSimulation(const EuroOption& option, const EuroOption& lower, const EuroOption& upper, double deltaS) {
    // GENERATE SAMPLES
    z_.resize(512 * 10000);
    
    std::generate(z_.begin(), z_.end(), BSM::standard_normal);
    
    auto zit = z_.cbegin();
    
    // PREPARE FOR MONTE CARLO
    std::array<double, 3> sum_of_outcome;
    std::fill(sum_of_outcome.begin(), sum_of_outcome.end(), 0.);
    
    // BEGINNING AND END OF SUBSEQUENCES
    std::vector<std::size_t> begins({0, 1, 2, 4, 8, 16, 32, 64, 128, 256});
    std::vector<std::size_t> ends({1, 2, 4, 8, 16, 32, 64, 128, 256, 512});
    
    std::cout << std::fixed << std::setprecision(6);
    
    for (std::size_t loop = 0; loop < begins.size(); loop++) {
        
        // GENERATE PRICES AND GREEKS FOR EVERY ITERATION
        for (std::size_t i = begins[loop] * 10000; i < ends[loop] * 10000; i++) {
            sum_of_outcome[0] += option.Call(*zit);
            sum_of_outcome[1] += lower.Call(*zit);
            sum_of_outcome[2] += upper.Call(*zit);
            zit++;
        }
    
        // FIND THE CURRENT MEAN (MONTE CARLO ESTIMATE)
        std::array<double, 3> mean;
        std::transform(sum_of_outcome.cbegin(), sum_of_outcome.cend(), mean.begin(), [&](double sum)->double {
            return sum / (ends[loop] * 10000.);
        });
        
        // OUTPUT RESULTS
        std::cout << "======= " << ends[loop] << " =======" << std::endl;
        std::cout << "Delta(C), error(Delta(C)), Gamma(C), error(Gamma(C))" << std::endl;
        
        double deltaC = (mean[2] - mean[1]) / (2. * deltaS);
        double gammaC = (mean[2] + mean[1] - mean[0] * 2.) / (deltaS * deltaS);
        std::cout << deltaC << std::endl;
        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(deltaC - option.DeltaCall()) << std::endl;
        std::cout << gammaC << std::endl;
        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(gammaC - option.GammaCall()) << std::endl;
    }
    
    
}
