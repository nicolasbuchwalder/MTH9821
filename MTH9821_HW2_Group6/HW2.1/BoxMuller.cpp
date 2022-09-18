// BoxMuller: static class that generates two independant normal (0,1) values from the Box Muller method, specifically the MarsagliañBray algorithm
// @ MTH9821 Homework2 Group6

#include "BoxMuller.hpp"

#include "LCG.hpp"

#include <cmath>

std::pair<double, double> BoxMuller::standard_normal_pair() {
    double u1, u2, X;

    do {
        u1 = 2. * LCG::gen() - 1.;
        u2 = 2. * LCG::gen() - 1.;
        X = u1 * u1 + u2 * u2;
    } while (X > 1);

    double Y = std::sqrt(-2. * std::log(X) / X);

    return std::make_pair(u1 * Y, u2 * Y);
}

//std::tuple<int, double, double> BoxMuller::standard_normal_pair_count() {
//    double u1, u2, X;
//    int uniform_count = 0;
//    do {
//        u1 = 2. * LCG::gen() - 1.;
//        u2 = 2. * LCG::gen() - 1.;
//        X = u1 * u1 + u2 * u2;
//        uniform_count += 2;
//    } while (X > 1);
//
//    double Y = std::sqrt(-2. * std::log(X) / X);
//
//    return std::make_tuple(uniform_count, u1 * Y, u2 * Y);
//}
