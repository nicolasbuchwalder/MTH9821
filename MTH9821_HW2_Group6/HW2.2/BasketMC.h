#pragma once

#include "EuropeanOption.h"

#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <numeric>



class BasketMC
{
private:

	double _S01;		// current price of the first underlying
	double _S02;		// current price of the second underlying
	double _sigma1;		// (constant) volatility of the first underlying
	double _sigma2;		// (constant) volatility of the second underlying
	double _r;		// continuous risk free rate (constant)
	double _T;		// time of the simulation
	double _rho;	// correlation between the two processes

	std::size_t _n;		// number of paths
	std::size_t _m;		// number of time intervals
	
	double _dt;		// time increments

	std::vector<double> _S1;		// array of the values of the first underlying
	std::vector<double> _S2;		// array of the values of the second underlying
	std::vector<double> _maxS1S2;	// maximum value of both underlying for each path

	std::vector<double> _C;		// array of the values of the contracts based on the underlying



public:
	// constructor
	BasketMC(std::size_t m, std::size_t n, double T, double S01, double S02, double sigma1, double sigma2, double r, double _rho);

	// main function of the class: runs all other subfunctions
	double launchSimulation(bool lookback, double K);

private:

	// calculates the next increment values for all paths
	void next();

	// functions to calculate the value for next time increment of a specific path  for both underlyings
	double nextS1(std::size_t i, double z1);
	double nextS2(std::size_t i, double z1, double z2);
	

	// pricing a path independant
	double PriceBasketIndependant(double K);

	// pricing a path dependant
	double PriceBasketDependant(double K);



};

