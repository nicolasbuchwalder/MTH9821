#include "EuropeanOption.h"
#include "BasketMC.h"
#include "BasketAnalysis.h"
#include "HestonMC.h"
#include "ImpliedVol.h"


int main() {

	std::cout << "========================================" << std::endl;
	std::cout << "=====  PART 1: VARIANCE REDUCTION ======" << std::endl;
	std::cout << "========================================" << std::endl << std::endl;



	std::cout << std::endl << std::endl;
	std::cout << "========================================" << std::endl;
	std::cout << "=== PART 2: PATH INDEP. BASKET OPTION ==" << std::endl;
	std::cout << "========================================" << std::endl << std::endl;

	BasketAnalysis ba1(1, 10000);

	for (int k = 0; k < 9; k++) {
		BasketMC mc2(1, ba1.new_iter(k), 0.5, 26, 29, 0.31, 0.21, 0.025, 0.3);
		ba1.print_iter(mc2.launchSimulation(false, 50));
	}
	
	std::cout << std::endl << std::endl;
	std::cout << "========================================" << std::endl;
	std::cout << "==== PART 3: PATH DEP. BASKET OPTION ===" << std::endl;
	std::cout << "========================================" << std::endl << std::endl;

	BasketAnalysis ba2(150, 50);

	for (int k = 0; k < 10; k++) {
		BasketMC mc2(150, ba2.new_iter(k), 0.5, 26, 29, 0.31, 0.21, 0.025, 0.3);
		ba2.print_iter(mc2.launchSimulation(true, 50));
	}

	std::cout << std::endl << std::endl;
	std::cout << "========================================" << std::endl;
	std::cout << "========= PART 4: HESTON MODEL =========" << std::endl;
	std::cout << "========================================" << std::endl << std::endl;

	EuropeanOption opt4(50, 50, 0.5, 0.3, 0.05, 0.);

	std::vector<std::size_t> num_paths4;
	std::vector<double> estimations4;

	for (int k = 0; k < 6; k++) {
		num_paths4.push_back(500 * std::pow(2, k));

		HestonMC mc4(num_paths4.back(), 175, 0.5, 50., 0.05, 0.09, 0.35 * 0.35, 4., 0.25, -0.15);

		estimations4.push_back(mc4.launchSimulation(opt4, false));
	}

	ImpliedVol impl(50, 50, 0.5, 0.05, 0., false, num_paths4, estimations4);
	impl.launch();

}