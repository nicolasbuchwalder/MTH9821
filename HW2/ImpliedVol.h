#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <tuple>
#include <numbers>

class ImpliedVol
{
private:

    double S0_;     // spot price
    double K_;      // strike price
    double T_;      // maturity
    double r_;      // interest rate (constant & coumpounded continuously)
    double q_;      // dividend rate (constant & coumpounded continuously)
    bool call_;     // boolean to say if option is a call or put

    double r_disc_; // discount by interest rate: exp(-r(T-t))
    double q_disc_; // discount by dividend rate: exp(-q(T-t))

    std::vector<std::size_t> num_paths;		// number of paths
	std::vector<double> prices;			// prices of the options
	std::vector<double> implied_vols;	// implied volatility of the prices

public:

    // constructor
    ImpliedVol(double S0_, double K_, double T_, double r_, double q_, bool call, std::vector<std::size_t> num_paths, std::vector<double> prices);

    void launch();


private:

    // computing the implied vol for ith price
    double Newton(std::size_t i, double sigma0 = 0.25, double tol = 0.000001);

    // printing the results
    void print();


    // price of an option in Black Scholes framework for a certain volatility
    double PriceBS(double sigma, double d1, double d2);

    // vega of an option in Black Scholes framework for a certain volatility
    double VegaBS(double sigma, double d1);

    // helper function to compute d1 and d2 in Black Scholes framework for a certain volatility
    std::pair<double, double> D1D2(double sigma);

    // helper function to compute cdf of normal
    double Phi(double z);


    

    

};

