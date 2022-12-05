//
//  LCG.hpp
//  MonteCarlo1
//
//  Created by Aubree Li on 11/22/22.
//

#ifndef LCG_h
#define LCG_h

class LCG {

private:
    // three parameters necessary for LCE
    const static thread_local unsigned long a_= 39373;
    const static thread_local unsigned long c_ = 0;
    const static thread_local unsigned long k_ = 2147483647;

    // starting state, can be seen as the seed of the generator
    static thread_local unsigned long state_;

public:
    // as we are only using static methods, we delete constructor
    LCG() = delete;
    ~LCG() = default;

    // static method that will generate values of pseudo uniform distribution
    static double gen() {
        
        // doing the mod seperately to reduce risk of overflow
        thread_local unsigned long state_ = 1;
        state_ = static_cast<unsigned long>(static_cast<unsigned long long>(a_) * static_cast<unsigned long long>(state_) % k_);
        state_ += c_ % k_;
        state_ %= k_;   // maybe ok to remove?

        return double(state_) / k_;
    }

    // static method to change the seed, which set the state parameter
    static void reseed(unsigned long new_seed) {
        state_ = new_seed;
    }
};

#endif /* LCG_h */
