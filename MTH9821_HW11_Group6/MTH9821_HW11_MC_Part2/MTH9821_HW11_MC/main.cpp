//
//  main.cpp
//  MonteCarlo1
//
//  Created by Aubree Li on 11/22/22.
//

#include <iostream>
#include "EuropeanOption.hpp"
#include "MonteCarlo.hpp"
#include <cmath>
#include <iomanip>

int main(int argc, const char * argv[]) {
    
    // problem 2 - No CV version
    EuropeanOption opt(50, 55.55, (7./12.), .2, .02, .0);
    
    // time increments list
    std::vector<double> dt_ls;
    dt_ls.push_back(2./12.);
    dt_ls.push_back(2./12.);
    dt_ls.push_back(2./12.);
    dt_ls.push_back(1./12.);
    
    std::cout << std::fixed << std::setprecision(6);
    
    // no cv, no variance reduction
    // number of paths simulated
    std::vector<int> N_no_cv;
    for (int k = 0; k<=8 ;++k) {
        N_no_cv.push_back(40000 * pow(2,k));
    }
    MonteCarlo mc_no_cv(opt, dt_ls, N_no_cv);
    mc_no_cv.Discrete_Div_MC(false);
    
    
    
    std::cout << "\nControl Variate:\n";
    // control variate
    std::vector<int> N_cv = N_no_cv;
    N_cv.pop_back();    // k =0:7 for cv
    // cv, variance reduction
    MonteCarlo mc_cv(opt, dt_ls, N_cv);
    mc_cv.Discrete_Div_MC(true);
    return 0;
}
