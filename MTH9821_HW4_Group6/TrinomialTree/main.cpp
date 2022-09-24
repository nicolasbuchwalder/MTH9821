//
//  main.cpp
//  BinomialTree
//
//  Created by Mingsen Wang on 2022/09/14.
//

#include <iostream>
#include "TrinomialTree.hpp"
#include "BTreeIV.hpp"
#include "BinomialTree.hpp"
#include "EuropeanOption.hpp"
#include <iomanip>

void PriceAndGreeks(const EuropeanOption& option, bool american) {
    std::cout << "\n\n\n" << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        TrinomialTree tree(41, .25, 1., i, .03, .005);
        auto res = tree.TreePricer(option.PutPayoff(), american);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }
    
    std::cout << "\n\n\n" << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        TrinomialTree tree(41, .25, 1., i, .03, .005);
        auto res = tree.TreePricerBS(option.PutPayoff(), american);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }
    
    std::cout << "\n\n\n" << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        TrinomialTree tree(41, .25, 1., i, .03, .005);
        auto res = tree.TreePricerBSR(option.PutPayoff(), american);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }
}

void PriceAndGreeks_compensated(const EuropeanOption& option) {
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        TrinomialTree tree(41, .25, 1., i, .03, .005);
        auto res = tree.TreePricer_PathDependent_compensated(option.PutPayoff(), option);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }
    
    std::cout << "\n\n\n" << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        TrinomialTree tree(41, .25, 1., i, .03, .005);
        auto res = tree.TreePricerBS_PathDependent_compensated(option.PutPayoff(), option);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }
    
    std::cout << "\n\n\n" << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        TrinomialTree tree(41, .25, 1., i, .03, .005);
        auto res = tree.TreePricerBSR_PathDependent_compensated(option.PutPayoff(), option);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }
}

void Problem12() {
    std::cout << std::setprecision(6);
    std::cout << std::fixed;
    std::cout << "======= Problem 1 ========" << std::endl;
    std::cout << "European - BS" << std::endl;
    EuropeanOption option(0., 41., 39., 1., .25, .03, .005);
    std::cout << option.Put() << '\t' << option.DeltaPut() << '\t' << std::setprecision(7) << option.GammaPut() << '\t' << std::setprecision(6) << option.ThetaPut() << std::endl;
    PriceAndGreeks(option, false);
    
    std::cout << "\n\n\n======= Problem 2 ========" << std::endl;
    std::cout << "American - Exact" << std::endl;
    TrinomialTree tree_exact1(41, .25, 1., 10000, .03, .005);
    TrinomialTree tree_exact2(41, .25, 1., 10001, .03, .005);
    
    auto exact1 = tree_exact1.TreePricer(option.PutPayoff(), true);
    auto exact2 = tree_exact2.TreePricer(option.PutPayoff(), true);
    
    std::cout << (exact1.value + exact2.value) / 2. << '\t' << (exact1.delta + exact2.delta) / 2. << '\t' << std::setprecision(7) << (exact1.gamma + exact2.gamma) / 2. << '\t' << std::setprecision(6) << (exact1.theta + exact2.theta) / 2. << std::endl;
    PriceAndGreeks(option, true);
    
    std::cout << "\n\n\n======= Variance Reduction =======" << std::endl;
    PriceAndGreeks_compensated(option);
}

void Problem3() {
    
    std::cout << std::setprecision(6);
    
    TTreeIV iv(48., 50., 1. / 3., .02, .01, 4.08);
    
    iv.PutIV();
}

int main(int argc, const char * argv[]) {
    std::cout << std::fixed << std::setprecision(6) << std::endl;

    Problem12();
    Problem3();
    return 0;
}
