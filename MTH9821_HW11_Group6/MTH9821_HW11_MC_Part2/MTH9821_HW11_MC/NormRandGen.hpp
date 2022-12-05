//
//  RandGen.hpp
//  MonteCarlo1
//
//  Created by Aubree Li on 11/22/22.
//

#ifndef NormRandGen_hpp
#define NormRandGen_hpp

#include <tuple>

// Inverse Transform Method
class ITM {

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
    // function that implements pseudo random normal distribution from BSM algorithm
    static double standard_normal();
};


// ---------------------------------------------
// Accept-Reject Method
class ARM
{
public:
    static double standard_normal();
    // generate standard normal while counting
    static std::tuple<int, double> standard_normal_count();
};


// ---------------------------------------------
// Box-Muller Method
class BoxMuller {
public:
    static std::pair<double, double> standard_normal_pair();
    static std::tuple<int, double, double> standard_normal_pair_count();
};

#endif /* NormRandGen_hpp */
