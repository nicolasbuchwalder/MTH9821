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


void Problem1_simulate(bool american) {
    EuropeanOption option(0., 54., 50., 1., .29, .0375, .01);
    
    for (std::size_t i = 10; i < 101; i+=1) {
        BinomialTree tree(54, .29, 1., i, .0375, .01);
        auto BT = tree.TreePricer(option.PutPayoff(), american);
        auto ABT = tree.AvgTreePricer(option.PutPayoff(), american);
        auto BBS = tree.TreePricerBS(option.PutPayoff(), american);
        auto BBSR = tree.TreePricerBSR(option.PutPayoff(), american);
        
        std::cout << BT.value << '\t' << ABT.value << '\t' << BBS.value << '\t' << BBSR.value << std::endl;
    }
}

void Problem1() {
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "======= Problem 1 ========" << std::endl;
    std::cout << "EUROPEAN" << std::endl;
    Problem1_simulate(false);
    std::cout << "AMERICAN" << std::endl;
    Problem1_simulate(true);
    std::cout << "\n\n";
}

void PriceAndGreeks(const EuropeanOption& option, bool american) {
    std::cout << "\n\n\n" << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        BinomialTree tree(54, .29, 1., i, .0375, .01);
        auto res = tree.TreePricer(option.PutPayoff(), american);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }
    
    std::cout << "\n\n\n" << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        BinomialTree tree(54, .29, 1., i, .0375, .01);
        auto res = tree.AvgTreePricer(option.PutPayoff(), american);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }
    
    std::cout << "\n\n\n" << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        BinomialTree tree(54, .29, 1., i, .0375, .01);
        auto res = tree.TreePricerBS(option.PutPayoff(), american);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }
    
    std::cout << "\n\n\n" << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        BinomialTree tree(54, .29, 1., i, .0375, .01);
        auto res = tree.TreePricerBSR(option.PutPayoff(), american);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }
}

void PriceAndGreeks_compensated(const EuropeanOption& option) {
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        BinomialTree tree(54, .29, 1., i, .0375, .01);
        auto res = tree.TreePricer_PathDependent_compensated(option.PutPayoff(), option);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }
    
    std::cout << "\n\n\n" << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        BinomialTree tree(54, .29, 1., i, .0375, .01);
        auto res = tree.AvgTreePricer_PathDependent_compensated(option.PutPayoff(), option);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }
    
    std::cout << "\n\n\n" << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        BinomialTree tree(54, .29, 1., i, .0375, .01);
        auto res = tree.TreePricerBS_PathDependent_compensated(option.PutPayoff(), option);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }
    
    std::cout << "\n\n\n" << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        BinomialTree tree(54, .29, 1., i, .0375, .01);
        auto res = tree.TreePricerBSR_PathDependent_compensated(option.PutPayoff(), option);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }
}

void Problem23() {
    std::cout << std::setprecision(6);
    std::cout << std::fixed;
    std::cout << "======= Problem 2 ========" << std::endl;
    std::cout << "European - BS" << std::endl;
    EuropeanOption option(0., 54., 50., 1., .29, .0375, .01);
    std::cout << option.Put() << '\t' << option.DeltaPut() << '\t' << std::setprecision(7) << option.GammaPut() << '\t' << std::setprecision(6) << option.ThetaPut() << std::endl;
    PriceAndGreeks(option, false);
    
    std::cout << "\n\n\n======= Problem 3 ========" << std::endl;
    std::cout << "American - Exact" << std::endl;
    BinomialTree tree(54, .29, 1., 10000, .0375, .01);
    auto exact = tree.AvgTreePricer(option.PutPayoff(), true);
    std::cout << exact.value << '\t' << exact.delta << '\t' << std::setprecision(7) << exact.gamma << '\t' << std::setprecision(6) << exact.theta << std::endl;
    PriceAndGreeks(option, true);
    
    std::cout << "\n\n\n======= Variance Reduction =======" << std::endl;
    PriceAndGreeks_compensated(option);
}

void Problem4() {
    
    std::cout << std::setprecision(6);
    
    BTreeIV iv(58., 60., .75, .02, .01, 6.36);
    
    iv.PutIV();
}

void Trinomial() {
    EuropeanOption option(0., 54., 50., 1., .29, .0375, .01);
    
    BinomialTree tree1(54, .29, 1., 400, .0375, .01);
    TrinomialTree tree2(54, .29, 1., 200, .0375, .01);
    auto res1 = tree2.TreePricer(option.PutPayoff(), true);
    auto res2 = tree2.TreePricer(option.PutPayoff(), true);
    
    std::cout << res1.value << std::endl;
    std::cout << res2.value << std::endl;
}

int main(int argc, const char * argv[]) {
    std::cout << std::fixed << std::setprecision(6) << std::endl;

    Problem1();
    
    Problem23();
    
    Problem4();

    // Trinomial();
    
    return 0;
}
