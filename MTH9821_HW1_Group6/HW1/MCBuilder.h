// MCBuilder: main class for Monte Carlo method estimation of an option
// taken out of the Advanced C++ course
// @ MTH9821 Homework1 Group6 

#ifndef MCBuilder_h
#define MCBuilder_h

#include "Mediator.h"

class MCBuilder {
private:
    double S0_;     // Initial price of the underlying asset
    double K_;      // Strike price
    double T_;      // Time to expiration
    double r_;      // Riskless rate
    double q_;      // Dividend rate
    double sigma_;  // Volatility (constant)
    
    unsigned NThreads_;         // # of threads
    unsigned long NSim_;        // # of simulations per thread
    unsigned NT_;               // # of time steps
    unsigned long Nk_;
    
    std::vector<std::shared_ptr<ISde>> sdes_;
    std::vector<std::shared_ptr<FdmBase>> fdms_;
    std::vector<std::shared_ptr<Mediator>> mediators_;
    std::vector<std::vector<std::shared_ptr<Pricer>>> pricers_;

    void CreateSde();
    void CreateFdm();
    void CreateMediator();
    void CreatePricer();

    double real_value;
    
public:
    // Get base option data
    MCBuilder(int k, bool optimal, const std::vector<double>& option, const double real_value);
    ~MCBuilder() = default;
    
    // Build the system in sequence
    void BuildMC();
    
    // Print the results
    void Print() const;
    
};

#endif /* MCBuilder_h */
