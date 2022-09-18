//
//  MultiStepMC.hpp
//  HW2.1
//
//  Created by 王明森 on 2022/09/10.
//

#ifndef MultiStepMC_hpp
#define MultiStepMC_hpp

#include <vector>
#include "LookbackBasketOption.hpp"

class MultiStepMC {
private:
    std::vector<std::vector<double>> z_;
    
public:
    MultiStepMC(std::size_t steps, std::size_t paths, std::size_t multiplier);
    ~MultiStepMC() = default;
    
    void Price(const LookbackBasketOption& option) const;
};

#endif /* MultiStepMC_hpp */
