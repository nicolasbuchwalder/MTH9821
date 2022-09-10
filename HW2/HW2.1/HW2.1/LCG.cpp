// LCG: static class that generates uniform (0,1) values from the Linear Congruential Generator method
// @ MTH9821 Homework2 Group6

#include "LCG.hpp"

// assigning parameters
thread_local unsigned long LCG::a_ = 39373;
thread_local unsigned long LCG::c_ = 0;
thread_local unsigned long LCG::k_ = 2147483647;    // 2^31-1
// assigning state
thread_local unsigned long LCG::state_ = 1;

// static method that will generate values of pseudo uniform distribution
double LCG::gen() {

    // doing the mod seperately to reduce risk of overflow
    state_ = static_cast<unsigned long>(static_cast<unsigned long long>(a_) * static_cast<unsigned long long>(state_) % k_);
    state_ += c_ % k_;
    state_ %= k_;

    return double(state_) / k_;
}

// static method to change the seed, which set the state parameter
void LCG::reseed(unsigned long new_seed) {
    state_ = new_seed;
}
