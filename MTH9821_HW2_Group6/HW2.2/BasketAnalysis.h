#pragma once

#include <vector>
#include <iostream>
#include <iomanip>

class BasketAnalysis
{
private:
	std::size_t _num_paths;	// number of paths of current iteration
	int _constant;			// constant that multiplies 2^k
	std::size_t _m;			// number of time intervals

public:

	// constructor
	BasketAnalysis(std::size_t m, int constant);

	
	size_t new_iter(size_t k);

	void print_iter(double price);
};

