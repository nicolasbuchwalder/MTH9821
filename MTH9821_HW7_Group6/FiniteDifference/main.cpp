//
//  main.cpp
//  MTH_9821_HW7_Group6
//
//  Created by 王明森 on 2022/10/21.
//

#include <iostream>
#include <iomanip>
#include "EuropeanOption.hpp"
#include "FiniteDifferencePricer.hpp"

void PrintVector(const std::vector<double>& vec) {
    for (auto elem : vec) {
        std::cout << elem << "\t";
    }
    std::cout << std::endl;
}

void Euro() {
    EuropeanOption option(0., 37., 40, .75, .28, .03, .015);
    FiniteDifferencePricer pricer(37., 40, .75, .28, .03, .015);
    
    std::cout << std::setprecision(8);
    std::cout << "\n\nEUROPEAN PUT FINITE DIFFERENCE" << std::endl;
    std::cout << "Exact: " << option.Put() << std::endl;

    for (std::size_t i = 4; i < 257; i <<= 2) {
        PrintVector(pricer.EuroPut(i));
    }
}

void American() {
    EuropeanOption option(0., 42., 40, .75, .32, .04, .02);
    FiniteDifferencePricer pricer(42., 40, .75, .32, .04, .02);
    
    std::cout << std::setprecision(8);
    std::cout << "\n\nAMERICAN PUT FINITE DIFFERENCE" << std::endl;
    std::cout << "Exact: " << 3.3045362802172642 << std::endl;
    
    for (std::size_t i = 4; i < 257; i <<= 2) {
        PrintVector(pricer.AmericanPut(i));
    }
}

void American_boundary() {
    FiniteDifferencePricer pricer(42., 40, .75, .32, .04, .02);
    
    std::cout << std::setprecision(8);
    std::cout << "\n\nAMERICAN PUT EARLY EX DOMAIN" << std::endl;
    
    std::vector<double> t;
    std::vector<double> Sopt;
    auto vecs = pricer.AmericanPut_EarlyExDomain(16);
    t = vecs[0];
    Sopt = vecs[1];
    PrintVector(t);
    PrintVector(Sopt);
}

int main(int argc, const char * argv[]) {
    Euro();
    American();
    American_boundary();
    
    return 0;
}
