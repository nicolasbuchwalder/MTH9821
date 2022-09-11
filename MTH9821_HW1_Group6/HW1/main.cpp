#include "EuropeanOption.h"
#include "OneStepMC.h"
#include "MCBuilder.h"

#include <iostream>
#include <vector>
#include <numeric>
#include <iomanip>
#include <cmath>
#include "LCG.h"
#include "ITM.h"
#include "ARM.h"
#include "BoxMuller.h"

// function that gives the price of a down and out option for ex3
double down_and_out() {
    EuropeanOption orig_option(39., 39., .75, .25, .02, .01);
    double B = 36.;
    EuropeanOption mod_option(B * B / 39., 39., .75, .25, .02, .01);

    double a = (.02 - .01) / (.25 * .25) - .5;

    return orig_option.Call() - pow((B / 39.), 2 * a) * mod_option.Call();
}




int main(int argc, const char * argv[]) {

    // PART 1
    // pricing path independant options
    std::cout << "========================================" << std::endl;
    std::cout << "==== PART 1:PATH INDEPENDANT OPTION ====" << std::endl;
    std::cout << "========================================" << std::endl << std::endl << std::endl;

    EuropeanOption option(41., 42., .75, .25, .03, .01);

    OneStepMC mc;

    mc.DirectSimulation(option);

    std::cout << std::endl;

    double deltaS = 0.01;
    EuropeanOption lower(41. - deltaS, 42., .75, .25, .03, .01);
    EuropeanOption upper(41. + deltaS, 42., .75, .25, .03, .01);
    LCG::reseed(1);
    mc.FDMSimulation(option, lower, upper, deltaS);

    std::cout << std::endl;


    // PART 2
    // pricing path dependant options

    std::cout << "========================================" << std::endl;
    std::cout << "===== PART 2:PATH DEPENDANT OPTION =====" << std::endl;
    std::cout << "========================================" << std::endl << std::endl << std::endl;

    std::vector<double> option2({ 39., 39., .75, .02, .01, .25 });

    double real_value = down_and_out();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Nk: \t\tm_k     n_k     C       err(C)"<< std::endl;
    for (int i = 0; i < 10; i++) {

        LCG::reseed(1);
        // fixed number of paths
        MCBuilder fixed_mont_carlo(i, false, option2, real_value);
        fixed_mont_carlo.BuildMC();
        fixed_mont_carlo.Print();

        LCG::reseed(1);
        // optimal number of paths and timestamps
        MCBuilder optimal_monte_carlo(i, true, option2, real_value);
        optimal_monte_carlo.BuildMC();
        optimal_monte_carlo.Print();

    }

    std::cout << std::endl;

    // PART 3
    // comparing esstimates for different random number generators

    std::cout << "========================================" << std::endl;
    std::cout << "==== PART3:DIFFERENT RANDOM ENGINES ====" << std::endl;
    std::cout << "========================================" << std::endl << std::endl << std::endl;

    EuropeanOption option3(50, 55., .5, .3, .04, .0);
    OneStepMC mc2;

    // ITM
    LCG::reseed(1);
    mc2.DiffRngSimulation(option3, 1);

    // ARM
    LCG::reseed(1);
    mc2.DiffRngSimulation(option3, 2);

    // BoxMuller
    LCG::reseed(1);
    mc2.DiffRngSimulation(option3, 3);



    return 0;
}
