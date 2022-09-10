// BoxMuller: static class that generates two independant normal (0,1) values from the Box Muller method, specifically the MarsagliañBray algorithm
// @ MTH9821 Homework2 Group6

#ifndef BOXMULLER_hpp
#define BOXMULLER_hpp

#include <tuple>

class BoxMuller {
public:
    static std::pair<double, double> standard_normal_pair();
    // static std::tuple<int, double, double> standard_normal_pair_count();
};
#endif /* BOXMULLER_hpp */
