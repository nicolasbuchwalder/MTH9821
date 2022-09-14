#pragma once

#include "EuropeanOption.h"

#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <numeric>



class HestonMC
{
private:

	double _S0;		// current price of the underlying
	double _mu;		// drift
	double _V0;		// initial variance of the underlying
	double _Vmean;	// long term variance mean	
	double _lambda; // speed of mean reversion 
	double _eta;	// standard deviation of the underlying variance
	double _rho;	// correlation between the two processes

	std::size_t _n;	// number of paths
	std::size_t _m; // number of time intervals
	double _T;		// time of the simulation
	double _dt;		// time increments

	std::vector<double> _V;		// array of the values of the varation of the returns
	std::vector<double> _S;		// array of the values of the underlying

	std::vector<double> _C;		// array of the values of the contracts based on the underlying
	std::vector<double> _Cimpl;	// array of the BlackScholes implied volatility of the contracts based on the underlying



public:
	// constructor
	HestonMC(std::size_t n, std::size_t m, double T, double S0, double mu, double V0, double Vmean, double lambda, double eta, double rho);

	// main function of the class: runs all other subfunctions
	double launchSimulation(EuropeanOption opt, bool call);

private:

	// calculates the next increment values for all paths
	void next();

	// functions to calculate the value for next time increment of a specific path (both variance and value)
	double nextV(std::size_t i, double z1, double z2);
	double nextS(std::size_t i, double z1);

	// pricing an european option
	double PriceEuropean(EuropeanOption opt, bool call);
	

};

