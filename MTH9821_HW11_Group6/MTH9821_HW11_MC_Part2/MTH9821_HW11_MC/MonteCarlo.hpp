//
//  MonteCarlo.hpp
//  MonteCarlo1
//
//  Created by Aubree Li on 11/22/22.
//

#ifndef MonteCarlo_hpp
#define MonteCarlo_hpp

#include <vector>
#include "EuropeanOption.hpp"

class MonteCarlo {
private:
    EuropeanOption opt;
    double _S0;        // current price of the first underlying
    double _sigma;        // (constant) volatility of the first underlying
    double _r;        // continuous risk free rate (constant)
    double _T;        // time of the simulation
    double _K;       // strike
    double _q;       // dividend rate
    
    std::vector<int> _N;        // list: number of paths
    std::size_t _m;        // number of time intervals
    std::vector<double> _dt;        // time increments

    std::vector<double> _S;        // array of the simulated S_T
    
public:
    MonteCarlo(const EuropeanOption& option, std::vector<double> dt_ls, std::vector<int> _N);
    
    // call this func for final result
    void Discrete_Div_MC(bool is_cv = false);
    
private:

    // functions to calculate the value for next time increment of a path for underlying asset
    double nextS(double S_prev, std::size_t i, double z, double prop_div=0, double fixed_div=0);
    

    // simulate max(_N) paths for S(T), stored in _S
    std::vector<double> Discrete_Div_simulation();

    
};

#endif /* MonteCarlo_hpp */
