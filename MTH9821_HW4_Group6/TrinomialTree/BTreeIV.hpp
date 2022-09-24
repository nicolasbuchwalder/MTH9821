//
//  BTreeIV.hpp
//  BinomialTree
//
//  Created by 王明森 on 2022/09/17.
//

#ifndef BTreeIV_hpp
#define BTreeIV_hpp

#include "BinomialTree.hpp"

class BTreeIV {
private:
    double S_;      // Spot price
    double K_;      // Strike price
    double T_;      // Maturity
    double r_;      // Const interest rate
    double q_;      // Dividend rate
    double V0_;     // Value at t=0
    
    std::size_t steps_; // Optimal steps
    
public:
    BTreeIV(double S, double K, double T, double r, double q, double V0);
    ~BTreeIV() = default;
    
    double PutIV() const;
};

class TTreeIV {
private:
    double S_;      // Spot price
    double K_;      // Strike price
    double T_;      // Maturity
    double r_;      // Const interest rate
    double q_;      // Dividend rate
    double V0_;     // Value at t=0
    
    std::size_t steps_; // Optimal steps
    
    static bool american_;  // American IV
    
public:
    TTreeIV(double S, double K, double T, double r, double q, double V0);
    ~TTreeIV() = default;
    
    double PutIV() const;
};

#endif /* BTreeIV_hpp */
