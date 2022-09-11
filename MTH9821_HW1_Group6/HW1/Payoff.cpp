// Payoff: utility classes used to calculate payoff of an option at maturity (put or call)
// taken out of the Advanced C++ course
// @ MTH9821 Homework1 Group6 


#include "Payoff.h"
#include <algorithm>

IPayoff::IPayoff(double K, bool call) : K_(K), Call_(call) {}

CallPayoff::CallPayoff(double K) : IPayoff(K, true) {}

double CallPayoff::operator () (double S) {
    return std::max(S - K_, 0.);
}

PutPayoff::PutPayoff(double K) : IPayoff(K, false) {}

double PutPayoff::operator () (double S) {
    return std::max(K_ - S, 0.);
}
