#include "ImpliedVol.h"

ImpliedVol::ImpliedVol(double S0, double K, double T, double r, double q, bool call, std::vector<std::size_t> num_paths, std::vector<double> prices)
	: S0_(S0), K_(K), T_(T), r_(r), q_(q), call_(call), num_paths(num_paths), prices(prices) 
{
	q_disc_ = std::exp(-q * T);
	r_disc_ = std::exp(-r * T);
};

void ImpliedVol::launch() {
	for (int i = 0; i < num_paths.size(); i++) {
		implied_vols.push_back(Newton(i));
	}
	print();
}


double ImpliedVol::Phi(double z) {
	return std::erfc(-z / std::sqrt(2.)) / 2.;
}

std::pair<double, double> ImpliedVol::D1D2(double sigma) {
	double d1 = (log(S0_ / K_) + (r_ - q_ + sigma * sigma / 2.) * T_) / (sigma * std::sqrt(T_));
	return std::make_pair(d1, d1 - sigma * std::sqrt(T_));
}


double ImpliedVol::PriceBS(double sigma, double d1, double d2) {
	if (call_) {
		return S0_ * q_disc_ * Phi(d1) - K_ * r_disc_ * Phi(d2);
	}
	else {
		return K_ * r_disc_ * Phi(-d2) - S0_ * q_disc_ * Phi(-d1);
	}
}

double ImpliedVol::VegaBS(double sigma, double d1) {
	return S0_* q_disc_* std::sqrt(T_) * std::exp(-d1 * d1 / 2.) / std::sqrt(2. * std::numbers::pi);
}

double ImpliedVol::Newton(std::size_t i, double sigma0, double tol) {
	std::pair<double, double> d1d2;

	double x_new = sigma0;
	double x_old = sigma0 + 10 * tol;

	while (std::abs(x_new - x_old) > tol) {
		x_old = x_new;
		d1d2 = D1D2(x_old);
		x_new = x_old - (PriceBS(x_old, d1d2.first, d1d2.second) - prices[i]) / VegaBS(x_old, d1d2.first);
	}
	return x_new;
}



void ImpliedVol::print() {
	std::cout << std::fixed << std::setprecision(6);
	std::cout << "N:    \tV       \tIV" << std::endl;
	for (int i = 0; i < num_paths.size(); i++) {
		std::cout << num_paths[i] << "\t";
		std::cout << prices[i] << ",  ";
		std::cout << implied_vols[i] << std::endl;
	}
}