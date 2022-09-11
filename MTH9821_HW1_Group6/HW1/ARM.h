// ARM: static class that generates normal (0,1) values from the Acceptance-Rejection method
// @ MTH9821 Homework1 Group6 

#ifndef ARM_h
#define ARM_h
#include <tuple>

class ARM
{
public:
    static double standard_normal();
    // generate standard normal while counting
    static std::tuple<int, double> standard_normal_count();
};

#endif /* ARM_h */