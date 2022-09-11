// MCBuilder: main class for Monte Carlo method estimation of an option
// taken out of the Advanced C++ course
// @ MTH9821 Homework1 Group6 


#include "MCBuilder.h"
#include <iostream>
#include <thread>
#include <numeric>

#include <cmath>


MCBuilder::MCBuilder(int k, bool optimal, const std::vector<double>& option, const double real_value) 
    : S0_(option[0]), K_(option[1]), T_(option[2]), r_(option[3]), q_(option[4]), sigma_(option[5]), Nk_(10000 * std::pow(2,k)), real_value(real_value) {
    
    // Get simulation parameters
//    std::cout << "# of threads? (1-8) ";
//    std::cin >> NThreads_;
//    NThreads_ %= 8; // Just in case someone tries to provide a weird number of threads
//    if (!NThreads_) NThreads_ += 8;
    NThreads_ = 1;

//    std::cout << "# of simulations? ";
//    std::cin >> NSim_;
//    NSim_ /= NThreads_; // Convert to simulations per thread
//    NSim_++;    // Just in case someone wants zero simulations
    

//    std::cout << "# of time steps? ";
//    std::cin >> NT_;
    if (optimal) {
        NT_ = std::ceil(pow(double(pow(2, k) * 10000.), 1. / 3) * pow(T_, 2. / 3));
        NSim_ = std::floor(pow(2, k) * 10000. / double(NT_));
    }
    else {
        NT_ = 200;
        NSim_ = 50 * pow(2, k);
    }
    
    // Build as many sde/fdm/mediators as threads
    sdes_ = std::vector<std::shared_ptr<ISde>>(NThreads_, nullptr);
    fdms_ = std::vector<std::shared_ptr<FdmBase>>(NThreads_, nullptr);
    mediators_ = std::vector<std::shared_ptr<Mediator>>(NThreads_, nullptr);
}

void MCBuilder::BuildMC() {
    
    // Get components sequentially
    this->CreateSde();
    this->CreateFdm();
    this->CreateMediator();
    
    // Allow multiple pricers for one simulation
    bool more_pricers = false;
//    bool more_pricers = true;
    do {
        // Add a pricer
        this->CreatePricer();

        // Ask if more pricers are needed
//        std::cout << "More pricers? 1: Yes 0: No " << std::endl;
//        std::cin >> more_pricers;

    } while (more_pricers);
    
    // Begin timing
//    std::cout << "\nRunning...\n" << std::endl;
    //StopWatch sw;
   // sw.StartStopWatch();

    // Simulate on different threads
    std::vector<std::thread> threads;
    std::for_each(mediators_.begin(), mediators_.end(), [&](auto& mediator) { threads.emplace_back(&Mediator::simulate, mediator); });

    // Finish up
    std::for_each(threads.begin(), threads.end(), [](auto& thread) { thread.join(); });
    //sw.StopStopWatch();
//    std::cout << "\nSimulation time: " << sw.GetTime() << " milliseconds.\n" << std::endl;
}

void MCBuilder::Print() const {
    for (const auto& identical_pricers : pricers_) {
        // Print type
//        identical_pricers.front()->PrintHeader();

        
        double estimated_price = std::accumulate(identical_pricers.cbegin(), identical_pricers.cend(), 0., [](double cum_price, const auto& pricer) {
            return cum_price + pricer->GetResults().price;
            }) / NThreads_;

        double estimated_variance = std::accumulate(identical_pricers.cbegin(), identical_pricers.cend(), 0., [](double cum_price, const auto& pricer) {
            return cum_price + pricer->GetResults().price;
            }) * (NSim_ - 1) / (NThreads_ * NSim_ - 1);

        double estimated_SD = std::sqrt(estimated_variance);
        double estimated_SE = std::sqrt(estimated_variance / NThreads_ / NSim_);
        std::cout << Nk_ << ":    \t";
        std::cout << NT_ << ", \t" << NSim_ << ", \t";
        std::cout << estimated_price << ", ";
        std::cout << std::abs(estimated_price - real_value) << ", ";
        std::cout << std::endl;
 

//        std::cout << "Estimated price   : " << estimated_price << std::endl;
//        std::cout << "Standard deviation: " << estimated_SD << std::endl;
//        std::cout << "Standard error    : " << estimated_SE << std::endl << std::endl;
    }
}

void MCBuilder::CreateSde() {
//    std::cout << "Choose SDE:" << std::endl;
//    std::cout << "0: Geometric Brownian motion" << std::endl;
//    std::cout << "1: Constant elasticity of variance" << std::endl;
    
    int SDEType = 0;
//    std::cin >> SDEType;
    
    double beta = 0;
    
    switch (SDEType) {
        case 0:
            // Geometric Brownian motion
            std::for_each(sdes_.begin(), sdes_.end(), [&](auto& x) { 
                x = std::make_shared<GBM>(r_, sigma_, q_, S0_, T_);
                }) ;
            break;
            
        case 1:
            // Constant elasticity of variance
            // We need one more parameter (beta)
            std::cout << "CEV parameter (beta): ";
            std::cin >> beta;
            std::for_each(sdes_.begin(), sdes_.end(), [&](std::shared_ptr<ISde>& x) { 
                x = std::make_shared<CEV>(r_, sigma_, q_, S0_, T_, beta); 
                });
            break;
            
        default:
            break;
    }
}

void MCBuilder::CreateFdm() {
//    std::cout << "Choose FDM:" << std::endl;
//    std::cout << "0: Explicit Euler" << std::endl;
//    std::cout << "1: Exact solution" << std::endl;
//    std::cout << "2: Milstein" << std::endl;
//    std::cout << "3: Predictor-corrector" << std::endl;
//    std::cout << "4: Modified predictor-corrector" << std::endl;
//    std::cout << "5: Multiplicative method" << std::endl;
    
    int FDMType = 5;
//    std::cin >> FDMType;
    
    // Parameters used in predictor-corrector methods
    double a = 0.;
    double b = 0.;

    switch (FDMType) {
        case 0:
            // Explicit Euler
            std::transform(sdes_.begin(), sdes_.end(), fdms_.begin(), [&](auto& sde) {
                return std::make_shared<Euler>(sde, NT_);
                });
            break;
            
        case 1:
            // Exact solution
            std::transform(sdes_.begin(), sdes_.end(), fdms_.begin(), [&](auto& sde) {
                return std::make_shared<Exact>(sde, NT_, sigma_, r_ - q_);
                });
            break;
            
        case 2:
            // Milstein method
            std::transform(sdes_.begin(), sdes_.end(), fdms_.begin(), [&](auto& sde) {
                return std::make_shared<Milstein>(sde, NT_);
                });
            break;
            
        case 3:
            // Predictor-corrector method
            std::cout << "Predictor-corrector parameters (separated by a whitespace):";
            std::cin >> a >> b;
            std::transform(sdes_.begin(), sdes_.end(), fdms_.begin(), [&](auto& sde) {
                return std::make_shared<PredictorCorrector>(sde, NT_, a, b);
                });
            break;
            
        case 4:
            // Modified predictor-corrector method
            std::cout << "Modified predictor-corrector parameters (separated by a whitespace):";
            std::cin >> a >> b;
            std::transform(sdes_.begin(), sdes_.end(), fdms_.begin(), [&](auto& sde) {
                return std::make_shared<ModifiedPredictorCorrector>(sde, NT_, a, b);
                });
            break;
            
        case 5:
            // Multiplicative method
            std::transform(sdes_.begin(), sdes_.end(), fdms_.begin(), [&](auto& sde) {
                return std::make_shared<Multiplicative>(sde, NT_, sigma_, r_ - q_);
                });
            break;
            
        default:
            break;
    }
}

void MCBuilder::CreateMediator() {
    auto sde_it = sdes_.begin();
    auto fdm_it = fdms_.begin();

    std::for_each(mediators_.begin(), mediators_.end(), [&](auto& x) {
        x = std::make_shared<Mediator>(*(sde_it++), *(fdm_it++), NSim_, NT_);
        });
}

void MCBuilder::CreatePricer() {
//    std::cout << "Choose option type:" << std::endl;
//    std::cout << "00: Vanilla European put" << std::endl;
//    std::cout << "01: Vanilla European call" << std::endl;
//    std::cout << "10: Asian put" << std::endl;
//    std::cout << "11: Asian call" << std::endl;
//    std::cout << "20: Barrier European put" << std::endl;
//    std::cout << "21: Barrier European call" << std::endl;
    
    int option_type = 21;
//    int option_type = 0;
//    std::cin >> option_type;
    
    // Used for barrier options
    double barrier = 0.;
    int barrier_type = 0;
    
    // Used for knock-out options
    double rebate = 0.;
    
    std::vector<std::shared_ptr<Pricer>> identical_pricers(NThreads_, nullptr);

    switch (option_type / 10) {
        case 0:
            // Vanilla European option
            // Note: Put/Call is determined by the last digit of option_type. Odd->Call; Even->Put
            //      The same trick is used for all other major option types.
            std::for_each(identical_pricers.begin(), identical_pricers.end(), [&](auto& x) {
                x = std::make_shared<EuropeanPricer>(option_type % 2, K_, std::exp(-1 * r_ * T_), NSim_);
                });
            break;
            
        case 1:
            // Asian option
            std::for_each(identical_pricers.begin(), identical_pricers.end(), [&](auto& x) {
                x = std::make_shared<AsianPricer>(option_type % 2, K_, std::exp(-1 * r_ * T_), NSim_);
                });
            break;
            
        case 2:
            // Barrier European call
            
            // Get barrier type
//            std::cout << "Select barrier option type:" << std::endl;
//            std::cout << "0: Down-and-out" << std::endl;
//            std::cout << "1: Down-and-in" << std::endl;
//            std::cout << "2: Up-and-out" << std::endl;
//            std::cout << "3: Up-and-in" << std::endl;
//            std::cin >> barrier_type;
            
            barrier_type = 0;
            // Note: If a number greater than 3 is provided, it is equivalent to providing the number mod 4.
            //      Similar to the trick that we determine if the option is a put/call.
            
            // Get the barrier
//            std::cout << "Barrier: ";
//            std::cin >> barrier;
            barrier = 36.;
            
            // If it is a knock-out option, ask for the rebate.
//            if (!(barrier_type % 2)) {
//                std::cout << "Rebate (if none, input 0): ";
//                std::cin >> rebate;
//            }
            
            std::for_each(identical_pricers.begin(), identical_pricers.end(), [&](auto& x) {
                x = std::make_shared<BarrierPricer>(option_type % 2, barrier_type / 2 % 2, barrier_type % 2, barrier, rebate, K_, std::exp(-1. * r_ * T_), NSim_);
                });
            break;
            
                
        default:
            break;
    }

    // Connect those parallel pricers to respective mediators
    auto mediator_it = mediators_.begin();
    std::for_each(identical_pricers.begin(), identical_pricers.end(), [&](auto& pricer) {
        (*mediator_it++)->AddPricer(pricer);
        });

    // Add this row of identical pricers to our system
    pricers_.emplace_back(std::move(identical_pricers));
}
