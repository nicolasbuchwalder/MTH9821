// OneStepMC: class that prices path independent options with Monte Carlo method
// @ MTH9821 Homework1 Group6 

#include "OneStepMC.h"
#include "EuropeanOption.h"
#include "ITM.h"
#include "ARM.h"
#include "BoxMuller.h"

#include <cmath>
#include <iostream>
#include <numbers>
#include <iomanip>
#include <algorithm>

// function to estimate price and greeks of an option with Monte Carlo
// with different random number generator
void OneStepMC::DirectSimulation(const EuropeanOption& option) {

    // GENERATE SAMPLES
    z_.resize(512 * 10000);

    std::generate(z_.begin(), z_.end(), ITM::standard_normal);

    auto zit = z_.cbegin();
   
    // PREPARE FOR MONTE CARLO
    std::array<double, 6> sum_of_outcome;
    std::fill(sum_of_outcome.begin(), sum_of_outcome.end(), 0.);

    // BEGINNING AND END OF SUBSEQUENCES
    std::vector<std::size_t> begins({ 0, 1, 2, 4, 8, 16, 32, 64, 128, 256 });
    std::vector<std::size_t> ends({ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512 });

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "======= DIRECT SIMULATION =======" << std::endl;
    std::cout << "N:\t\tC         error(C)  Delta(C)  error(Delta(C)) \tVega(C)    err(Vega(C))\tP         error(P)   Delta(P)  err(Delta(P))\tVega(P)    err(Vega(P))" << std::endl;

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
        
        std::cout << ends[loop] * 10000 <<":     \t";

        // call results
        std::cout << mean[0] << ", ";
        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[0] - option.Call()) << ", ";
        std::cout << mean[1] << ", ";
        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[1] - option.DeltaCall()) << ",  \t";
        std::cout << mean[2] << ", ";
        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[2] - option.VegaCall()) << ",  \t";

        // put results
        std::cout << mean[3] << ", ";
        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[3] - option.Put()) << ", ";
        std::cout << mean[4] << ", ";
        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[4] - option.DeltaPut()) << ", \t";
        std::cout << mean[5] << ", ";
        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(mean[5] - option.VegaPut()) << std::endl;
    }

}


// function to estimate price of an option with Monte Carlo and greeks with finite differences method
// with different random number generator
void OneStepMC::FDMSimulation(const EuropeanOption& option, 
    const EuropeanOption& lower, 
    const EuropeanOption& upper, 
    double deltaS) {
    // GENERATE SAMPLES
    z_.resize(512 * 10000);

    std::generate(z_.begin(), z_.end(), ITM::standard_normal);
    std::cout << "=======   FSD SIMULATION   =======" << std::endl;
    std::cout << "N:\t\tDelta(C)  error(Delta(C))  Gamma(C)  err(Gamma(C))"<< std::endl;
    auto zit = z_.cbegin();

    // PREPARE FOR MONTE CARLO
    std::array<double, 3> sum_of_outcome;
    std::fill(sum_of_outcome.begin(), sum_of_outcome.end(), 0.);

    // BEGINNING AND END OF SUBSEQUENCES
    std::vector<std::size_t> begins({ 0, 1, 2, 4, 8, 16, 32, 64, 128, 256 });
    std::vector<std::size_t> ends({ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512 });

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
        std::cout << ends[loop] * 10000 << ":     \t";

        double deltaC = (mean[2] - mean[1]) / (2. * deltaS);
        double gammaC = (mean[2] + mean[1] - mean[0] * 2.) / (deltaS * deltaS);
        std::cout << deltaC << ", ";
        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(deltaC - option.DeltaCall()) << ",\t   ";
        std::cout << gammaC << ", ";
        std::cout << std::sqrt(10000. * ends[loop]) * std::abs(gammaC - option.GammaCall()) << std::endl;
    }


}


// Part 3 compare different rngs
void OneStepMC::DiffRngSimulation(const EuropeanOption& option, int generator_id)
{
    // generator_id choices: 1: "Inverse Transform", 2: "Accept-Rejection", 3: "Box-Muller"
    
    // GENERATE SAMPLES with different rngs
    z_.resize(512 * 10000, 0);
    std::vector<double> M;      // record the number of independent std normal samples generated given diff uniform count
    std::vector<int> N{ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512 };
    int uniform_count;

    // generate standard normal samples using different rngs
    switch (generator_id)
    {
    default:
        break;
    case 1:
        // use Inverse Transform rng
        std::generate(z_.begin(), z_.end(), ITM::standard_normal);
        // M independent std normal samples
        M = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512 };
        std::cout << "======= INVERSE TRANSFORM METHOD =======" << std::endl;
        break;
    case 2:
        uniform_count = 0;
        // use accept-rejection rng
        z_.resize(0);
        for (auto& count_thresh : N)
        {
            do
            {
                auto arm_rv = ARM::standard_normal_count();
                z_.push_back(std::get<1>(arm_rv));
                
                uniform_count += std::get<0>(arm_rv);
            } while (uniform_count < (count_thresh * 10000.));

            // when thresh is hitted, summarize the number of standard normal rv generated
            M.push_back(double(z_.size())/10000);
        }
        std::cout << "==== ACCEPTANCE REJRECTION METHOD ====" << std::endl;
        break;

    case 3:
        // use Box-Muller
        uniform_count = 0;
        z_.resize(0);
        for (auto& count_thresh : N)
        {
            do
            {
                auto bm_rv = BoxMuller::standard_normal_pair_count();
                z_.push_back( std::get<1>(bm_rv));
                z_.push_back(std::get<2>(bm_rv));
                uniform_count += std::get<0>(bm_rv);
            } while (uniform_count < count_thresh * 10000);

            // when thresh is hitted, summarize the number of standard normal rv generated
            M.push_back(double(z_.size()) / 10000);
        }
        std::cout << "======= BOX MUELLER METHOD =======" << std::endl;
        break;

    };

    auto zit = z_.cbegin();
    double sum_of_outcome = 0;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "N:\t\tcount   \tC         error(C)" << std::endl;

    int i_path = 0;

    for (std::size_t loop = 0; loop < M.size(); loop++) {
        
        // GENERATE PRICES FOR EVERY ITERATION
        while (i_path < (M[loop] * 10000.)) {
            std::array<double, 2> out;
            sum_of_outcome += option.Put(*zit);
            zit++;
            i_path++;
        }

        // FIND THE CURRENT MEAN (MONTE CARLO ESTIMATE)
        double mean = sum_of_outcome / (i_path+1);

        // OUTPUT RESULTS
       
        std::cout << static_cast<int>(N[loop] * 10000) << ":     \t";
        std::cout << static_cast<int>(M[loop] * 10000) << ",    \t";
        std::cout << mean << ", ";
        std::cout << std::abs(mean - option.Put()) << std::endl;

    }
    std::cout << std::endl;
}
