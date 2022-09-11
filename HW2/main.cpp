#include "EuropeanOption.h"
#include "HestonMC.h"
#include "ImpliedVol.h"

int main() {

	EuropeanOption opt3(50, 50, 0.5, 0.3, 0.05, 0.);

	std::vector<std::size_t> num_paths;
	std::vector<double> estimations;

	for (int k = 0; k < 6; k++) {
		num_paths.push_back(500 * std::pow(2, k));

		HestonMC mc3(num_paths.back(), 175, 0.5, 50., 0.05, 0.09, std::sqrt(0.35), 4., 0.25, -0.15);

		estimations.push_back(mc3.launchSimulation(opt3, false));
	}

	ImpliedVol impl(50, 50, 0.5, 0.05, 0., false, num_paths, estimations);
	impl.launch();

}