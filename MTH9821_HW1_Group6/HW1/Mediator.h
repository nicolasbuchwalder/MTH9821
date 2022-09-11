//
//  Mediator.hpp
//  Final_Project
//
//  Created by Wanderers' Library on 2022/2/26.
//

#ifndef Mediator_h
#define Mediator_h

#include "ITM.h"
#include <boost/signals2/signal.hpp>
#include <vector>
#include <memory>
#include <iostream>
#include "Sde.h"
#include "Fdm.h"
#include "Pricer.h"

class MIS {
public:
    void operator() (unsigned long i) const {
        std::cout << i << " iterations" << std::endl;
    }
};

class Mediator {
private:
    std::shared_ptr<ISde> sde_;
    std::shared_ptr<FdmBase> fdm_;
    // Uses a static RNG
    
    unsigned long NSim_;        // # of simulations
    unsigned NT_;               // # of time steps
    
    
    
    boost::signals2::signal<void(unsigned long)> progression_signal_;   // Signal the MIS after a certain number of loops
    boost::signals2::signal<void (std::vector<double>)> pricer_signal_; // Send path to pricer
    boost::signals2::signal<void ()> stop_signal_;                      // Tell pricer to calculate results
    boost::signals2::signal<void ()> print_signal_; // Tell pricer to print results
    
public:
    Mediator(const std::shared_ptr<ISde>& sde, const std::shared_ptr<FdmBase>& fdm, unsigned long NSim, unsigned NT);
    ~Mediator() = default;
    
    // Add observers
    void AddMis(const MIS& mis);            // Add an observing MIS
    void AddPricer(const std::shared_ptr<Pricer>& pricer);   // Add an observing pricer
    
    // Main method for path generation
    void simulate();
    
    void PrintResults() const;
};

#endif /* Mediator_h */
