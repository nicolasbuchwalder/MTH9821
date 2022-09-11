#include "HestonMC.h"

#include "BoxMuller.h"

HestonMC::HestonMC(std::size_t n, std::size_t m, double T, double S0, double mu, double V0, double Vmean, double lambda, double eta, double rho)
	: _n{ n }, _m{ m }, _T{ T }, _S0{ S0 }, _mu{ mu }, _V0{ V0 }, _Vmean{ Vmean }, _lambda{ lambda }, _eta{ eta }, _rho{ rho }
{
	_dt = _T / m;
	for (int i = 0; i < _n; i++) {
		_S.push_back(_S0);
		_V.push_back(_V0);
	}

}

double HestonMC::nextS(std::size_t i, double z1) {
	return _S[i] * std::exp((_mu - _V[i] / 2.) * _dt + std::sqrt(_V[i] * _dt) * z1);
}

double HestonMC::nextV(std::size_t i, double z1, double z2) {
	return std::max(_V[i] - _lambda * (_V[i] - _Vmean) * _dt + _eta * std::sqrt(_V[i] * _dt) * (_rho * z1 + std::sqrt(1 - _rho * _rho)) * z2, 0.);
};


void HestonMC::next() {
	for (int i = 0; i < _n; i++) {
		std::pair<double, double> z = BoxMuller::standard_normal_pair();
		_S[i] = nextS(i, z.first);
		_V[i] = nextV(i, z.first, z.second);
	}
}

double HestonMC::PriceEuropean(EuropeanOption opt, bool call) {
	if (call) {
		for (auto S : _S) {
			_C.push_back(opt.CallFromStock(S));
		}
	}
	else {
		for (auto S : _S) {
			_C.push_back(opt.PutFromStock(S));
		}
	}
	return std::accumulate(_C.begin(), _C.end(), 0.) / _C.size();
}

double HestonMC::launchSimulation(EuropeanOption opt, bool call) {

	for (int i = 0; i < _m; i++) {
		HestonMC::next();
	}

	return PriceEuropean(opt, call);

}