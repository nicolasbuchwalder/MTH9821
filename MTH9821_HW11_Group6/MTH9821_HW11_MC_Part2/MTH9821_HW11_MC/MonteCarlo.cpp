//
//  MonteCarlo.cpp
//  MonteCarlo1
//
//  Created by Aubree Li on 11/22/22.
//

#include "MonteCarlo.hpp"
#include "NormRandGen.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>

MonteCarlo::MonteCarlo(const EuropeanOption& option, std::vector<double> dt_ls, std::vector<int> n) : opt(option), _dt(dt_ls), _N(n), _S0(option.S0_), _sigma(option.sigma_), _r(option.r_), _T(option.T_), _K(option.K_), _q(option.q_) {}


double MonteCarlo::nextS(double S_prev, std::size_t i, double z, double prop_div, double fixed_div) {
    // add discrete dividend to estimate path
    return S_prev * std::exp((_r - _sigma * _sigma / 2.) * _dt[i] + _sigma * std::sqrt(_dt[i]) * z) * (1 - prop_div) - fixed_div;
}


std::vector<double> MonteCarlo::Discrete_Div_simulation () {
    // simulate max(_N) paths for S(T), stored in _S
    // plus non-dividend paying simulation
    std::vector<double> S_no_div;
    
    for (int n=0; n<_N[_N.size()-1]; n++) {
        std::vector<int> path_i;
        double S_temp, S_no_div_temp;
        
        // generate 4 z_i for the simulation, i.e. two paris with BoxMuller
        std::pair<double, double> z1 = BoxMuller::standard_normal_pair();
        std::pair<double, double> z2 = BoxMuller::standard_normal_pair();
        // t1: 75 cents fixed div
        S_temp = _S0;
        S_temp = nextS(S_temp, 0, z1.first, 0., 0.75);
        // t2: 0.02 proportional div
        S_temp = nextS(S_temp, 1, z1.second, 0.02, 0.);
        // t3: 25 cents fixed div
        S_temp = nextS(S_temp, 2, z2.first, 0., 0.25);
        // t4: no div, get S_T
        S_temp = nextS(S_temp, 3, z2.first, 0., 0.);
        
        // add S_temp to _S
        _S.push_back(S_temp);
        
        // control variate simulation
        // i.e. non-dividend paying asset simulation
        S_no_div_temp = _S0 * std::exp((_r - _sigma * _sigma / 2.)* _T  +
        _sigma * (std::sqrt(_dt[0]) * z1.first +
                  std::sqrt(_dt[1]) * z1.second +
                  std::sqrt(_dt[2]) * z2.first +
                  std::sqrt(_dt[3]) * z2.second));
        S_no_div.push_back(S_no_div_temp);
    }
    
    return S_no_div;
}
 

void MonteCarlo::Discrete_Div_MC(bool is_cv) {
    // simulate for S_T and control variate S_no_div
    std::vector<double> S_no_div = Discrete_Div_simulation();
    
    std::vector<double> P_div, P_no_div;
    // P_i, simulated with discrete div
    std::transform(_S.begin(), _S.end(), std::back_inserter(P_div), [=](double s) { return std::exp(-1.*_r*_T) * std::max((_K - s), 0.); });
    // P_no_dic_i
    std::transform(S_no_div.begin(), S_no_div.end(), std::back_inserter(P_no_div), [this](double s) { return std::exp(-1*_r*_T) * std::max((_K - s), 0.); });
    
    // delta_i, with dividend
    // std::vector<double> delta_div = _S;
    // std::transform(_S.begin(), _S.end(), delta_div.begin(), [this](double s) { return -1*std::exp(-1*_r*_T) * std::max((_K - s), 0.)  /(_K - s) * (s / _S0) * (1 - 0.01); });
    std::vector<double> delta_div;
    // stick to the original version from Dan's hw
    for (int i = 0; i < _S.size(); ++i) {
        double d_temp = -1 * std::exp(-1 * _r * _T) * std::max((_K - _S[i]), 0.) / (_K - _S[i]) * (S_no_div[i] / _S0) * (1 - 0.01);
        delta_div.push_back(d_temp);
    }

    // delta_no_div_i, with no dividend, cv
    std::vector<double> delta_no_div;
    std::transform(S_no_div.begin(), S_no_div.end(), std::back_inserter(delta_no_div), [this](double s) { return -1 * std::exp(-1*_r*_T) * std::max((_K - s), 0.) / (_K - s) * (s / _S0); });
    
    
    
    // no control variate version
    if (!is_cv) {
        // estimate under differen number of paths simulated
        for (auto& n : _N) {
            double P = std::accumulate(P_div.begin(), P_div.begin() + n, 0.) / n;
            
            double delta = std::accumulate(delta_div.begin(), delta_div.begin() + n, 0.) / n;
            
            // print result
            std::cout << n << "\t" << P << "\t" << delta << "\n";
        }
        
    } else {
        // control variate version
        double cov_sum_p = 0, var_sum_p = 0;
        double cov_sum_delta = 0, var_sum_delta = 0;
        
        for (int i = 0; i < _N.size(); ++i) {
            int n = _N[i];
            
            // mean for P
            double P = std::accumulate(P_div.begin(), P_div.begin() + n, 0.) / n;
            double P_cv = std::accumulate(P_no_div.begin(), P_no_div.begin() + n, 0.) / n;
            
            // mean for delta
            double delta = std::accumulate(delta_div.begin(), delta_div.begin() + n, 0.) / n;
            double delta_cv = std::accumulate(delta_no_div.begin(), delta_no_div.begin() + n, 0.) / n;
            
            // cumsum based on previous n, calculate the corr coefficient
            double start = (i == 0) ? 0 : _N[i-1];
            // calculate coefficient b for P and delta
            for (int j = start; j < n; j++) {
                cov_sum_p +=  (P_div[j] - P) * (P_no_div[j] - P_cv);
                var_sum_p += pow((P_no_div[j] - P_cv), 2);
                
                cov_sum_delta = (delta_div[j] - delta) * (delta_no_div[j] - delta_cv);
                var_sum_delta = pow((delta_no_div[j] - delta_cv), 2);
            }
            
            double b_p = cov_sum_p / var_sum_p;
            double b_delta = cov_sum_delta / var_sum_delta;

            // W_i for p and delta
            double P_BS = opt.Put();
            double delta_BS = opt.DeltaPut();
            
            double w_p_sum = 0, w_delta_sum = 0;
            for (int j=0 ;j<n; ++j) {
                w_p_sum += P_div[j] - b_p * (P_no_div[j] - P_BS);
                w_delta_sum += delta_div[j] - b_delta * (delta_no_div[j] - delta_BS);
            }
            
            std::cout << n << "\t" << (w_p_sum / (n * 1.)) << "\t" << (w_delta_sum / (n * 1.)) << "\n";
        }
    }
    
}
