// ITM: static class that generates normal (0,1) values from the Inverse Transform method, specifically the Beasley-Springer-Moro algorithm 
// @ MTH9821 Homework1 Group6 

#include "ITM.h"

#include "LCG.h"

#include <cmath>

// define constants
thread_local const double ITM::a0_ = 2.50662823884;
thread_local const double ITM::a1_ = -18.61500062529;
thread_local const double ITM::a2_ = 41.39119773534;
thread_local const double ITM::a3_ = -25.44106049637;

thread_local const double ITM::b0_ = -8.47351093090;
thread_local const double ITM::b1_ = 23.08336743743;
thread_local const double ITM::b2_ = -21.06224101826;
thread_local const double ITM::b3_ = 3.13082909833;

thread_local const double ITM::c0_ = 0.3374754822726147;
thread_local const double ITM::c1_ = 0.9761690190917186;
thread_local const double ITM::c2_ = 0.1607979714918209;
thread_local const double ITM::c3_ = 0.0276438810333863;
thread_local const double ITM::c4_ = 0.0038405729373609;
thread_local const double ITM::c5_ = 0.0003951896511919;
thread_local const double ITM::c6_ = 0.0000321767881768;
thread_local const double ITM::c7_ = 0.0000002888167364;
thread_local const double ITM::c8_ = 0.0000003960315187;

// function that implements pseudo random normal distribution from BSM algorithm
double ITM::standard_normal() {

    double u = LCG::gen();

    double y = u - 0.5;
    double x;
    if (std::abs(y) < 0.42) {
        double r = y * y;
        x = y * (((a3_ * r + a2_) * r + a1_) * r + a0_) / ((((b3_ * r + b2_) * r + b1_) * r + b0_) * r + 1);
    }
    else {
        double r = u;
        if (y > 0) r = 1. - u;

        r = std::log(-std::log(r));
        x = c0_ + r * (c1_ + r * (c2_ + r * (c3_ + r * (c4_ + r * (c5_ + r * (c6_ + r * (c7_ + r * c8_)))))));

        if (y < 0) x = -x;
    }
    return x;
}