// Pricer: class to calulate price the options through Monte Carlo methods
// taken out of the Advanced C++ course
// @ MTH9821 Homework1 Group6 

#ifndef Pricer_h
#define Pricer_h

#include <vector>
#include <memory>
#include "Payoff.h"

struct MCResults {
    double price;   // Estimated price
    double variance;// Sample variance
};

// Abstract one-factor option pricer
class Pricer {
protected:
    std::shared_ptr<IPayoff> payoff_;   // Call/Put
    double discounter_;                 // Discounting factor
    unsigned long NSim_;                // # of simulations
    std::vector<double> prices_;        // Simulated prices
    MCResults results_;                 // Result: price, SE
    
public:
    Pricer(bool is_call, double strike, double discounter, unsigned long NSim);
    virtual ~Pricer() = default;
    
    // Note: since pricers and the mediator run on different threads,
    //      pricers accept copies of paths.
    virtual void ProcessPath(std::vector<double> path) = 0;
    virtual void PostProcess() final;           // No need to be overridden
    virtual MCResults GetResults() const final; // No need to be overridden
    
    // Print header (e.g. "Vanilla European call")
    virtual void PrintHeader() const = 0;

    // Print price, SD, SE
    virtual void PrintResults() const = 0;
    
};

// Vanilla European option pricer
class EuropeanPricer : public Pricer {
public:
    EuropeanPricer(bool is_call, double strike, double discounter, unsigned long NSim);
    virtual ~EuropeanPricer() override final = default;
    
    virtual void ProcessPath(std::vector<double> path) override final;
    
    // Print option type, price, SD, SE
    virtual void PrintHeader() const override final;
    virtual void PrintResults() const override final;
};

// Asian option pricer
// Explainer: The "settlement" price is the average price of the underlying asset (the meaning of "average" is specified in the contract)
// This implementation implements "the geometric average of the latest 10% of the price path".
class AsianPricer : public Pricer {
private:
    std::function<double (const std::vector<double>&)> average_;
    
public:
    // Default average implementation
    AsianPricer(bool is_call, double strike, double discounter, unsigned long NSim);
    // User-provided average implementation
    AsianPricer(bool is_call, double strike, double discounter, unsigned long NSim, const std::function<double (const std::vector<double>&)>& average);
    virtual ~AsianPricer() override final = default;
    
    virtual void ProcessPath(std::vector<double> path) override final;
    
    virtual void PrintHeader() const override final;
    virtual void PrintResults() const override final;
};

// Barrier European option pricer
class BarrierPricer : public Pricer {
private:
    double barrier_;
    bool is_up_;
    bool is_in_;
    double rebate_;
public:
    // Two additional params: up/down and in/out
    BarrierPricer(bool is_call, bool is_up, bool is_in, double barrier, double strike, double discounter, unsigned NSim);
    
    // With knock-out rebate
    BarrierPricer(bool is_call, bool is_up, bool is_in, double barrier, double rebate, double strike, double discounter, unsigned NSim);
    virtual ~BarrierPricer() override final = default;
    
    virtual void ProcessPath(std::vector<double> path) override final;
    
    virtual void PrintHeader() const override final;
    virtual void PrintResults() const override final;
};



#endif /* Pricer_h */
