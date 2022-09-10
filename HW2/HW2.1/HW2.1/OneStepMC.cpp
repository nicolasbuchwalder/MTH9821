// OneStepMC: class that prices path independent options with Monte Carlo method
// @ MTH9821 Homework1 Group6

#include "OneStepMC.hpp"
#include "BoxMuller.hpp"
#include "LCG.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>

OneStepMC::OneStepMC(std::size_t size) {
    // Generate samples using Box Muller
    LCG::reseed(1); // Reset LCE
    z_.resize(size);
    auto zgen_it = z_.begin();
    
    while (zgen_it != z_.end()) {
        double z1, z2;
        std::tie(z1, z2) = BoxMuller::standard_normal_pair();
        *(zgen_it++) = z1;
        *(zgen_it++) = z2;
    }
}

void OneStepMC::ControlVariate(const EuropeanOption& option) const {
    // Find the BS value of the option
    double BS = option.Put();
    
    // Prepare for Monte Carlo
    double S_sum = 0.;
    double S2_sum = 0.;
    double V_sum = 0.;
    double SV_sum = 0.;
    
    // BEGINNING AND END OF SUBSEQUENCES
    std::vector<std::size_t> begins({ 0, 1, 2, 4, 8, 16, 32, 64, 128, 256 });
    std::vector<std::size_t> ends({ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512 });
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "======= CONTROL VARIATE =======" << std::endl;
//    std::cout << "N:\t\tC         error(C)  Delta(C)  error(Delta(C)) \tVega(C)    err(Vega(C))\tP         error(P)   Delta(P)  err(Delta(P))\tVega(P)    err(Vega(P))" << std::endl;
    
    auto zit = z_.cbegin();
    
    for (std::size_t loop = 0; loop < begins.size(); loop++) {

        // GENERATE PRICES AND GREEKS FOR EVERY ITERATION
        for (std::size_t i = begins[loop] * 10000; i < ends[loop] * 10000; i++) {
            double S_i, V_i;

            std::tie(S_i, V_i) = option.SpotAndPut(*(zit++));
            
            // Sum of outcome: S, S^2, V, S*V
            S_sum += S_i;
            S2_sum += S_i * S_i;
            V_sum += V_i;
            SV_sum += S_i * V_i;
        }

        double N = ends[loop] * 10000;
        
        // FIND THE CURRENT MEAN (MONTE CARLO ESTIMATE)
        double b_hat = (SV_sum - S_sum * V_sum / N) / (S2_sum - S_sum * S_sum / N);
        
        double W_hat = (V_sum / N) - b_hat * (S_sum / N - option.S0_ / option.r_disc_);

        // OUTPUT RESULTS

        std::cout << std::setw(7) << ends[loop] * 10000 <<": ";

        // call results
        std::cout << std::setw(8) << W_hat << ", " << std::setw(8) << std::abs(W_hat - BS) << std::endl;
    }
    
    std::cout << std::endl;
}

void OneStepMC::AntitheticVariable(const EuropeanOption& option) const {
    // Find the BS value of the option
    double BS = option.Put();
    
    // Prepare for Monte Carlo
    double V_sum = 0.;
    
    // BEGINNING AND END OF SUBSEQUENCES
    std::vector<std::size_t> begins({ 0, 1, 2, 4, 8, 16, 32, 64, 128, 256 });
    std::vector<std::size_t> ends({ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512 });
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "======= ANTITHETIC VARIABLE =======" << std::endl;
    
    auto zit = z_.cbegin();
    
    for (std::size_t loop = 0; loop < begins.size(); loop++) {

        // GENERATE PRICES AND GREEKS FOR EVERY ITERATION
        for (std::size_t i = begins[loop] * 10000; i < ends[loop] * 10000; i++) {
            
            double V1 = option.SpotAndPut(*zit)[1];
            double V2 = option.SpotAndPut(-(*zit))[1];
            
            V_sum += (V1 + V2) / 2.;
            
            zit++;
        }

        double N = ends[loop] * 10000;
        
        // FIND THE CURRENT MEAN (MONTE CARLO ESTIMATE)
        double V_hat = V_sum / N;

        // OUTPUT RESULTS

        std::cout << std::setw(7) << ends[loop] * 10000 <<": ";

        // call results
        std::cout << std::setw(8) << V_hat << ", " << std::setw(8) << std::abs(V_hat - BS) << std::endl;
    }
    
    std::cout << std::endl;
}

void OneStepMC::MomentMatching(const EuropeanOption& option) const {
    // Find the BS value of the option
    double BS = option.Put();
    
    // Prepare for Monte Carlo
    std::vector<double> S(z_.size());
    
    // Generate array of S
    auto GenerateS = [&](double z)->double {
        return option.SpotAndPut(z)[0];
    };
    std::transform(z_.cbegin(), z_.cend(), S.begin(), GenerateS);
    
    // BEGINNING AND END OF SUBSEQUENCES
    std::vector<std::size_t> begins({ 0, 1, 2, 4, 8, 16, 32, 64, 128, 256 });
    std::vector<std::size_t> ends({ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512 });
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "======= MOMENT MATCHING =======" << std::endl;
    
    double S_sum = 0.;
    
    for (std::size_t loop = 0; loop < begins.size(); loop++) {

        S_sum += std::accumulate((S.cbegin() + begins[loop] * 10000), (S.cbegin() + ends[loop] * 10000), 0.);
        double moment_factor = option.S0_ / option.r_disc_ / S_sum * (ends[loop] * 10000.);
        
        double V_sum = 0.;
        
        for (auto Sit = S.cbegin(); Sit != (S.cbegin() + ends[loop] * 10000); Sit++) {
            V_sum += (option.r_disc_ * std::max(option.K_ - (*Sit) * moment_factor, 0.));
        }
        
        // FIND THE CURRENT MEAN (MONTE CARLO ESTIMATE)
        double V_hat = V_sum / (ends[loop] * 10000.);

        // OUTPUT RESULTS

        std::cout << std::setw(7) << ends[loop] * 10000 <<": ";

        // call results
        std::cout << std::setw(8) << V_hat << ", " << std::setw(8) << std::abs(V_hat - BS) << std::endl;
    }
    
    std::cout << std::endl;
}

void OneStepMC::MMCV(const EuropeanOption& option) const {
    // Find the BS value of the option
    double BS = option.Put();
    
    // Prepare for Monte Carlo
    std::vector<double> S(z_.size());
    
    // Generate array of S
    auto GenerateS = [&](double z)->double {
        return option.SpotAndPut(z)[0];
    };
    std::transform(z_.cbegin(), z_.cend(), S.begin(), GenerateS);
    
    // BEGINNING AND END OF SUBSEQUENCES
    std::vector<std::size_t> begins({ 0, 1, 2, 4, 8, 16, 32, 64, 128, 256 });
    std::vector<std::size_t> ends({ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512 });
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "= MOMENT MATCHING and CONTROL VARIATES =" << std::endl;
    
    double S_sum = 0.;
    
    for (std::size_t loop = 0; loop < begins.size(); loop++) {

        S_sum += std::accumulate((S.cbegin() + begins[loop] * 10000), (S.cbegin() + ends[loop] * 10000), 0.);
        double moment_factor = option.S0_ / option.r_disc_ / S_sum * (ends[loop] * 10000.);
        
        // Generate moment-matched V
        std::vector<double> tildeV(ends[loop] * 10000);
        auto GenerateTildeV = [&](double S) -> double {
            return option.r_disc_ * std::max(option.K_ - S * moment_factor, 0.);
        };
        std::transform(S.cbegin(), (S.cbegin() + ends[loop] * 10000), tildeV.begin(), GenerateTildeV);
        
        double tildeStildeV_sum = 0.;
        double tildeV_sum = 0.;
        double tildeS2_sum = 0.;
        double tildeS_sum = S_sum * moment_factor;
        // Calculate sums for b_hat and W_hat
        auto Sit = S.cbegin();
        auto tildeVit = tildeV.cbegin();
        for (std::size_t i = 0; i < ends[loop] * 10000; i++) {
            
            tildeStildeV_sum += (*Sit) * (*tildeVit);
            tildeV_sum += *tildeVit;
            tildeS2_sum += (*Sit) * (*Sit);
            
            Sit++;
            tildeVit++;
        }
        
        tildeStildeV_sum *= moment_factor;
        tildeS2_sum *= (moment_factor) * (moment_factor);

        double N = ends[loop] * 10000;
        
        // FIND THE CURRENT MEAN (MONTE CARLO ESTIMATE)
        double b_hat = (tildeStildeV_sum - tildeV_sum * tildeS_sum) / (tildeS2_sum - 2. * tildeS_sum * option.S0_ / option.r_disc_ + N * option.S0_ / option.r_disc_);
        
        double W_hat = (tildeV_sum / N) - b_hat * (tildeS_sum / N - option.S0_ / option.r_disc_);
        
        // OUTPUT RESULTS
        std::cout << std::setw(7) << ends[loop] * 10000 <<": ";

        // call results
        std::cout << std::setw(8) << W_hat << ", " << std::setw(8) << std::abs(W_hat - BS) << std::endl;
    }
    
    std::cout << std::endl;
}

//// function to estimate price and greeks of an option with Monte Carlo
//// with different random number generator
//void OneStepMC::DirectSimulation(const EuropeanOption& option) {
//
//    // GENERATE SAMPLES
//    z_.resize(512 * 10000);
//
//    std::generate(z_.begin(), z_.end(), ITM::standard_normal);
//
//    auto zit = z_.cbegin();
//
//    // PREPARE FOR MONTE CARLO
//    std::array<double, 6> sum_of_outcome;
//    std::fill(sum_of_outcome.begin(), sum_of_outcome.end(), 0.);
//
//    // BEGINNING AND END OF SUBSEQUENCES
//    std::vector<std::size_t> begins({ 0, 1, 2, 4, 8, 16, 32, 64, 128, 256 });
//    std::vector<std::size_t> ends({ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512 });
//
//    std::cout << std::fixed << std::setprecision(6);
//    std::cout << "======= DIRECT SIMULATION =======" << std::endl;
//    std::cout << "N:\t\tC         error(C)  Delta(C)  error(Delta(C)) \tVega(C)    err(Vega(C))\tP         error(P)   Delta(P)  err(Delta(P))\tVega(P)    err(Vega(P))" << std::endl;
//
//    for (std::size_t loop = 0; loop < begins.size(); loop++) {
//
//        // GENERATE PRICES AND GREEKS FOR EVERY ITERATION
//        for (std::size_t i = begins[loop] * 10000; i < ends[loop] * 10000; i++) {
//            std::array<double, 6> out = option.Outcome(*(zit++));
//
//            for (std::size_t j = 0; j < 6; j++) {
//                sum_of_outcome[j] += out[j];
//            }
//
//        }
//
//        // FIND THE CURRENT MEAN (MONTE CARLO ESTIMATE)
//        std::array<double, 6> mean;
//        std::transform(sum_of_outcome.cbegin(), sum_of_outcome.cend(), mean.begin(), [&](double sum)->double {
//            return sum / (ends[loop] * 10000.);
//            });
//
//        // OUTPUT RESULTS
//
//        std::cout << ends[loop] * 10000 <<":     \t";
//
//        // call results
//        std::cout << mean[0] << ", ";
//        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[0] - option.Call()) << ", ";
//        std::cout << mean[1] << ", ";
//        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[1] - option.DeltaCall()) << ",  \t";
//        std::cout << mean[2] << ", ";
//        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[2] - option.VegaCall()) << ",  \t";
//
//        // put results
//        std::cout << mean[3] << ", ";
//        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[3] - option.Put()) << ", ";
//        std::cout << mean[4] << ", ";
//        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[4] - option.DeltaPut()) << ", \t";
//        std::cout << mean[5] << ", ";
//        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[5] - option.VegaPut()) << std::endl;
//    }
//
//}
//
