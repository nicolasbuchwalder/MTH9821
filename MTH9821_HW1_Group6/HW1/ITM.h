// ITM: static class that generates normal (0,1) values from the Inverse Transform method, specifically the Beasley-Springer-Moro algorithm 
// @ MTH9821 Homework1 Group6 
#ifndef ITM_hpp
#define ITM_hpp

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

#endif /*ITM_hpp*/ 
