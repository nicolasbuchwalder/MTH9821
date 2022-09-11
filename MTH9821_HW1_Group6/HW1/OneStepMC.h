// OneStepMC: class that prices path independent options with Monte Carlo method
// @ MTH9821 Homework1 Group6 
//
// Last Modifications
// 2022.09.05 add DiffRngSimulation

#ifndef OneStepMC_hpp
#define OneStepMC_hpp

#include <vector>
#include <string>

// forward definition of the class
class EuropeanOption;

class OneStepMC {
private:
    std::vector<double> z_;

public:

    OneStepMC() = default;
    ~OneStepMC() = default;

    void DirectSimulation(const EuropeanOption& option);
    void FDMSimulation(const EuropeanOption& option, const EuropeanOption& lower, const EuropeanOption& upper, double deltaS);
    void DiffRngSimulation(const EuropeanOption& option, int generator_id);     // generator_id choices: 1: "Inverse Transform", 2: "Accept-Rejection", 3: "Box-Muller"

};

#endif /* OneStepMC_hpp */
