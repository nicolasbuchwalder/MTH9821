//
//  PseudoRNG.hpp
//  HW1
//
//  Created by 王明森 on 2022/08/30.
//

#ifndef PseudoRNG_hpp
#define PseudoRNG_hpp

#include <utility>
#include <random>

class LCE_uniform {
    // std::linear_congruential_engine + std::uniform_real_distribution(0., 1.)
    // Custom implementation
private:
    static thread_local unsigned long a_;
    static thread_local unsigned long c_;
    static thread_local unsigned long k_;
    
    static thread_local unsigned long state_;
    
    static thread_local std::random_device rv_;
    static thread_local std::mt19937 eng_;
    static thread_local std::uniform_real_distribution<double> unif_dist_;
    
public:
    // No instances of LCE_uniform is needed.
    LCE_uniform() = delete;
    ~LCE_uniform() = default;
    
    // Generate new number
    static double gen();
    
    static void reseed(unsigned long new_seed);
};

class BSM {
    // Beasley-Springer-Moro algorithm
private:
    static thread_local const double a0_;
    static thread_local const double a1_;
    static thread_local const double a2_;
    static thread_local const double a3_;
    
    static thread_local const double b0_;
    static thread_local const double b1_;
    static thread_local const double b2_;
    static thread_local const double b3_;
    
    static thread_local const double c0_;
    static thread_local const double c1_;
    static thread_local const double c2_;
    static thread_local const double c3_;
    static thread_local const double c4_;
    static thread_local const double c5_;
    static thread_local const double c6_;
    static thread_local const double c7_;
    static thread_local const double c8_;
    
public:
    static double standard_normal();
};

class Rejection {
public:
    static double standard_normal();
};

class MarsagliaBray {
public:
    static std::pair<double, double> standard_normal_pair();
};

class Rng{
    // Standard normal distribution generator using Mersenne Twister
    // For testing and debugging only
    
private:
    // Random device
    static thread_local std::random_device rv_;
    
    // Random engine initially seeded with rv
    static thread_local std::mt19937_64 eng_;
    
    // Standard normal distribution
    static thread_local std::normal_distribution<double> normal_;
    
public:
    // No instances of Rng is needed.
    Rng() = delete;
    ~Rng() = default;
    
    // Generate new standard normal variable
    static double gen();
    
    static void reseed(unsigned new_seed);
};

#endif /* PseudoRNG_hpp */
