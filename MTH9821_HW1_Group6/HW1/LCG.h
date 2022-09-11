// LCG: static class that generates uniform (0,1) values from the Linear Congruential Generator method
// @ MTH9821 Homework1 Group6 

#ifndef LCG_hpp
#define LCG_hpp

class LCG {

private:
    // three parameters necessary for LCE
    static thread_local unsigned long a_;
    static thread_local unsigned long c_;
    static thread_local unsigned long k_;

    // starting state, can be seen as the seed of the generator
    static thread_local unsigned long state_;

public:
    // as we are only using static methods, we delete constructor
    LCG() = delete;
    ~LCG() = default;

    // static method that will generate values of pseudo uniform distribution
    static double gen();

    // static method to change the seed, which set the state parameter
    static void reseed(unsigned long new_seed);
};

#endif /* LCG_hpp */

