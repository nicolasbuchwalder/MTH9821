// OneStepMC: class that prices path independent options with Monte Carlo method
// @ MTH9821 Homework2 Group6
//
// Last Modifications
// 2022.09.05 add DiffRngSimulation

#ifndef OneStepMC_hpp
#define OneStepMC_hpp

#include <vector>
#include "EuropeanOption.hpp"

// forward definition of the class
class EuropeanOption;

class OneStepMC {
private:
    std::vector<double> z_;

public:
    OneStepMC(std::size_t size);
    ~OneStepMC() = default;
    
    void ControlVariate(const EuropeanOption& option) const;
    void AntitheticVariable(const EuropeanOption& option) const;
    void MomentMatching(const EuropeanOption& option) const;
    void MMCV(const EuropeanOption& option) const;

};

#endif /* OneStepMC_hpp */
