// Payoff: utility classes used to calculate payoff of an option at maturity (put or call)
// taken out of the Advanced C++ course
// @ MTH9821 Homework1 Group6 


#ifndef Payoff_h
#define Payoff_h

/*
 Utility classes used in pricer:
 Call payoff (max{S - K, 0})
 Put  payoff (max{K - S, 0})
 */

#include <functional>

class IPayoff {
protected:
    double K_;  // Strike price
    
public:
    const bool Call_; // 1 = Call, 0 = Put
    
    IPayoff(double K, bool call);
    virtual ~IPayoff() = default;
    
    virtual double operator() (double S) = 0;
};

class CallPayoff : public IPayoff {
public:
    CallPayoff(double K);
    virtual ~CallPayoff() override final = default;
    
    virtual double operator() (double S) override final;
};

class PutPayoff : public IPayoff {
public:
    PutPayoff(double K);
    virtual ~PutPayoff() override final = default;
    
    virtual double operator() (double S) override final;
};



#endif /* Payoff_h */
