#include "BasketAnalysis.h"

BasketAnalysis::BasketAnalysis(std::size_t m, int constant)
	: _m{ m }, _constant{ constant }, _num_paths{0}
{
	std::cout << std::fixed << std::setprecision(6);
	std::cout << "N:    \tm\tV" << std::endl;
};

size_t BasketAnalysis::new_iter(size_t k) {
	_num_paths = _constant * std::pow(2, k);
	return _num_paths;
}
void BasketAnalysis::print_iter(double price) {
	std::cout << _num_paths << "\t";	
	std::cout << _m << ",\t";
	std::cout << price << std::endl;
}