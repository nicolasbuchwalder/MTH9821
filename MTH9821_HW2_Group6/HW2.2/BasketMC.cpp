#include "BasketMC.h"

#include "BoxMuller.h"

BasketMC::BasketMC(std::size_t m, std::size_t n, double T, double S01, double S02, double sigma1, double sigma2, double r, double rho)
	: _m{ m }, _n{ n }, _T{ T }, _S01{ S01 }, _S02{ S02 }, _sigma1{ sigma1 }, _sigma2{ sigma2 }, _r{ r }, _rho{ rho }
{
	LCG::reseed(1);
	_dt = _T / _m;

	for (int i = 0; i < _n; i++) {
		_S1.push_back(_S01);
		_S2.push_back(_S02);

		_maxS1S2.push_back(_S01 + _S02);

	}
}
	
double BasketMC::nextS1(std::size_t i, double z1) {
	return _S1[i] * std::exp((_r - _sigma1 * _sigma1 / 2.) * _dt + _sigma1 * std::sqrt(_dt) * z1);
}

double BasketMC::nextS2(std::size_t i, double z1, double z2) {
	return  _S2[i] * std::exp((_r - _sigma2 * _sigma2 / 2.) * _dt + _sigma2 * std::sqrt(_dt) * (_rho * z1 + std::sqrt(1 - _rho * _rho) * z2));
};


void BasketMC::next() {
	for (int i = 0; i < _n; i++) {
		std::pair<double, double> z = BoxMuller::standard_normal_pair();

		_S1[i] = nextS1(i, z.first);
		_S2[i] = nextS2(i, z.first, z.second);

		_maxS1S2[i] = std::max(_maxS1S2[i], _S1[i] + _S2[i]);
	}
}

double BasketMC::PriceBasketIndependant(double K) {
	for (int i = 0; i < _n; i++) {
		_C.push_back(std::max(_S1[i] + _S2[i] - K, 0.));
	}
	return std::exp(-_r * _T) * std::accumulate(_C.begin(), _C.end(), 0.) / _C.size();
}

double BasketMC::PriceBasketDependant(double K) {
	for (int i = 0; i < _n; i++) {
		_C.push_back(std::max(_maxS1S2[i] - K, 0.));
	}
	return std::exp(-_r * _T) * std::accumulate(_C.begin(), _C.end(), 0.) / _C.size();
}

double BasketMC::launchSimulation(bool lookback, double K) {

	for (int i = 0; i < _m; i++) {
		BasketMC::next();
	}

	if (!lookback) {
		return PriceBasketIndependant(K);
	}
	else {
		return PriceBasketDependant(K);
	}
}