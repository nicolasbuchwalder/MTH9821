//
//  PseudoRNG.cpp
//  HW1
//
//  Created by 王明森 on 2022/08/30.
//

#include "PseudoRNG.hpp"
#include <cmath>


thread_local std::random_device LCE_uniform::rv_ = std::random_device();
thread_local std::mt19937 LCE_uniform::eng_ = std::mt19937(LCE_uniform::rv_());
thread_local std::uniform_real_distribution<double> LCE_uniform::unif_dist_ = std::uniform_real_distribution<double>(0., 1.);

// class LCE_uniform
// Parameters: x <- (a * x + c) (mod k)
thread_local unsigned long LCE_uniform::a_ = 39373;
thread_local unsigned long LCE_uniform::c_ = 0;
thread_local unsigned long LCE_uniform::k_ = 2147483647;    // 2^31-1

// Initial state
thread_local unsigned long LCE_uniform::state_ = 1;

double LCE_uniform::gen() {
    
//    return (LCE_uniform::unif_dist_)(LCE_uniform::eng_);
    
    state_ = a_ * state_ % k_;
    state_ += c_ % k_;
    state_ %= k_;
    
    return double(state_) / k_;
}

void LCE_uniform::reseed(unsigned long new_seed) {
    state_ = new_seed;
}

thread_local const double BSM::a0_ =   2.50662823884;
thread_local const double BSM::a1_ = -18.61500062529;
thread_local const double BSM::a2_ =  41.39119773534;
thread_local const double BSM::a3_ = -25.44106049637;

thread_local const double BSM::b0_ =  -8.47351093090;
thread_local const double BSM::b1_ =  23.08336743743;
thread_local const double BSM::b2_ = -21.06224101826;
thread_local const double BSM::b3_ =   3.13082909833;

thread_local const double BSM::c0_ = 0.3374754822726147;
thread_local const double BSM::c1_ = 0.9761690190917186;
thread_local const double BSM::c2_ = 0.1607979714918209;
thread_local const double BSM::c3_ = 0.0276438810333863;
thread_local const double BSM::c4_ = 0.0038405729373609;
thread_local const double BSM::c5_ = 0.0003951896511919;
thread_local const double BSM::c6_ = 0.0000321767881768;
thread_local const double BSM::c7_ = 0.0000002888167364;
thread_local const double BSM::c8_ = 0.0000003960315187;

double BSM::standard_normal() {
    
    double u = LCE_uniform::gen();
    
    double y = u - 0.5;
    double x;
    if (std::abs(y) < 0.42) {
        double r = y * y;
        x = y * (((a3_ * r + a2_) * r + a1_) * r + a0_) / ((((b3_ * r + b2_) * r + b1_) * r + b0_) * r + 1);
    } else {
        double r = u;
        if (y > 0) r = 1. - u;
        
        r = std::log(-std::log(r));
        x = c0_ + r * (c1_ + r * (c2_ + r * (c3_ + r * (c4_ + r * (c5_ + r * (c6_ + r * (c7_ + r * c8_)))))));
        
        if (y < 0) x = -x;
    }
    return x;
}

double Rejection::standard_normal() {
    // Generate three uniform rvs on [0, 1]
    double u1 = LCE_uniform::gen();
    double u2 = LCE_uniform::gen();
    double u3 = LCE_uniform::gen();
    
    // Generate an exponential rv with param 1
    double x = -std::log(u1);
    
    if (u2 > std::exp(-(x - 1.) * (x - 1.) / 2.)) {
        // Reject
        return Rejection::standard_normal();
    } else {
        // Accept
        if (u3 < .5) x = -x;    // Flip sign with .5 prob
        return x;
    }
}

std::pair<double, double> MarsagliaBray::standard_normal_pair() {
    double u1, u2, X;
    
    do {
        u1 = 2. * LCE_uniform::gen() - 1.;
        u2 = 2. * LCE_uniform::gen() - 1.;
        X = u1 * u1 + u2 * u2;
    } while (X > 1);
    
    double Y = std::sqrt(-2. * std::log(X) / X);
    
    return std::make_pair(u1 * Y, u2 * Y);
}

// Random device (used to seed the engine)
thread_local std::random_device Rng::rv_ = std::random_device();

// Seed engine with std::random_device
thread_local std::mt19937_64 Rng::eng_ = std::mt19937_64(Rng::rv_());

// Standard normal distribution N(0, 1)
thread_local std::normal_distribution<double> Rng::normal_ = std::normal_distribution<double>(0., 1.);

// Generate a new standard normal variable
double Rng::gen() {
    return ((Rng::normal_)(Rng::eng_));
}

// Re-seed the engine with custom seed
void Rng::reseed(unsigned new_seed) {
    Rng::eng_.seed(new_seed);
}
