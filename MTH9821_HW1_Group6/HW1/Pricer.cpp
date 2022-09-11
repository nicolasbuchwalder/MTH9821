// Pricer: class to calulate price the options through Monte Carlo methods
// taken out of the Advanced C++ course
// @ MTH9821 Homework1 Group6 

#include "Pricer.h"
#include <numeric>
#include <iostream>

/* Abstract one-factor option pricer */

Pricer::Pricer(bool is_call, double strike, double discounter, unsigned long NSim) : discounter_(discounter), NSim_(NSim) {
    
    if (is_call) {
        payoff_ = std::make_shared<CallPayoff>(strike);
    } else {
        payoff_ = std::make_shared<PutPayoff>(strike);
    }
    
    // There will be as many simulated prices as simulations
    prices_.reserve(NSim);
}

void Pricer::PostProcess() {
    // Note: Discounters are applied to results only.
    // So discounting occurs only twice instead of NSim times.
    
    double N = NSim_;
    
    results_.price = std::accumulate(prices_.cbegin(), prices_.cend(), 0.) / N;    // Expected price = average simulated price
    
    // SD^2 = (sum of squares - square of sum) / (N - 1)
    double sum_of_squares = std::accumulate(prices_.cbegin(), prices_.cend(), 0., [](double init, double x)->double {
        return init + x * x;
    });
    
    results_.variance = (sum_of_squares - results_.price * results_.price * N) / (N - 1);
    
    // Discount prices and variance!
    results_.price *= discounter_;
    results_.variance *= (discounter_ * discounter_);
}

MCResults Pricer::GetResults() const {
    return results_;
}

void Pricer::PrintResults() const {
    std::cout << "Estimated price: " << results_.price << std::endl;
    std::cout << "Standard dev.:   " << std::sqrt(results_.variance) << std::endl;
    std::cout << "Standard err.:   " << std::sqrt(results_.variance / NSim_) << std::endl;
}

/* Vanilla European option pricer */
EuropeanPricer::EuropeanPricer(bool is_call, double strike, double discounter, unsigned long NSim) : Pricer(is_call, strike, discounter, NSim) {}

void EuropeanPricer::ProcessPath(std::vector<double> path) {
    // Payoff is path-independent.
    // Undiscounted price = PayoffFunction(Asset price at expiry)
    prices_.push_back(payoff_->operator()(path.back()));
}

void EuropeanPricer::PrintHeader() const {
    // "Vanilla European call" or "Vanilla European put"
    std::cout << "Vanilla European ";
    std::cout << (payoff_->Call_ ? "call" : "put") << std::endl;
}

void EuropeanPricer::PrintResults() const {
    this->PrintHeader();
    
    // Print price, SD, SE
    this->Pricer::PrintResults();
}

/* Asian option pricer */
// This c'tor uses the default "average" implementation.
AsianPricer::AsianPricer(bool is_call, double strike, double discounter, unsigned long NSim) : Pricer(is_call, strike, discounter, NSim), average_([](const std::vector<double>& path) -> double {
    
    // "The geometric average of the latest 10% of the price path"
    long price_start = 9 * path.size() / 10;
    double product = std::accumulate(path.cbegin() + price_start, path.cend(), 1., std::multiplies<>());
    return std::pow(product, 1. / (path.size() - price_start));
    
}) {}

// Pricer with user-given average function
// Not actually used in the final implementation though
AsianPricer::AsianPricer(bool is_call, double strike, double discounter, unsigned long NSim, const std::function<double (const std::vector<double>&)>& average) : Pricer(is_call, strike, discounter, NSim), average_(average) {}

void AsianPricer::ProcessPath(std::vector<double> path) {
    // Settlement price
    double settlement = average_(path);
    prices_.push_back(payoff_->operator()(settlement));
}

void AsianPricer::PrintHeader() const {
    // "Asian call" or "Asian put"
    std::cout << "Asian ";
    std::cout << (payoff_->Call_ ? "call" : "put") << std::endl;
}

void AsianPricer::PrintResults() const {
    this->PrintHeader();
    
    // Print price, SD, SE
    this->Pricer::PrintResults();
}

/* Barrier European option pricer */
BarrierPricer::BarrierPricer(bool is_call, bool is_up, bool is_in, double barrier, double strike, double discounter, unsigned NSim) : Pricer(is_call, strike, discounter, NSim), is_up_(is_up), is_in_(is_in), barrier_(barrier), rebate_(0.) {}

BarrierPricer::BarrierPricer(bool is_call, bool is_up, bool is_in, double barrier, double rebate, double strike, double discounter, unsigned NSim) : Pricer(is_call, strike, discounter, NSim), is_up_(is_up), is_in_(is_in), barrier_(barrier), rebate_(rebate) {}

void BarrierPricer::ProcessPath(std::vector<double> path) {
    
    bool crossed = false;
    
    if (is_up_) {
        // Find a node below the barrier
        // Unfortunately, std::lower_bound cannot be used, as it uses binary search and requires the container to be sorted.
        auto below_barrier = std::find_if(path.cbegin(), path.cend(), [&](double price) {
            return price < barrier_;
            // Or we can use std::less with std::bind, but that would be too unintelligible.
        });
        
        // An up-crossing happens only if there is a node above the barrier after a node below the barrier
        auto over_barrier = std::find_if(below_barrier, path.cend(), [&](double price) {
            return price > barrier_;
            // Or we can use std::more with std::bind, but that would be too unintelligible.
        });
        
        crossed = (over_barrier != path.cend());
        
    } else {
        // Find a node over the barrier
        // Unfortunately, std::lower_bound cannot be used, as it uses binary search and requires the container to be sorted.
        auto over_barrier = std::find_if(path.cbegin(), path.cend(), [&](double price) {
            return price > barrier_;
            // Or we can use std::less with std::bind, but that would be too unintelligible.
        });
        
        // An down-crossing happens only if there is a node below the barrier after a node over the barrier
        auto below_barrier = std::find_if(over_barrier, path.cend(), [&](double price) {
            return price < barrier_;
            // Or we can use std::more with std::bind, but that would be too unintelligible.
        });
        
        crossed = (below_barrier != path.cend());
    }
    
    
    
    if (is_in_) {
        if (crossed) {
            // and-in and barrier crossed
            prices_.push_back(payoff_->operator()(path.back()));
        } else {
            // and-in but barrier not crossed
            prices_.push_back(0.);
        }
    } else {
        if (crossed) {
            // and-out but barrier crossed
            prices_.push_back(rebate_);
        } else {
            // and-out and barrier not crossed
            prices_.push_back(payoff_->operator()(path.back()));
        }
    }
}

void BarrierPricer::PrintHeader() const {
    // "Up/Down-and-in/out European call/put"
    std::cout << (is_up_ ? "Up" : "Down") << "-and-";
    std::cout << (is_in_ ? "in" : "out") << " European ";
    std::cout << (payoff_->Call_ ? "call" : "put") << std::endl;
}

void BarrierPricer::PrintResults() const {
    this->PrintHeader();
    
    // Print price, SD, SE
    this->Pricer::PrintResults();
}
