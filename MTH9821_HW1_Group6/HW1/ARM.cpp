// ARM: static class that generates normal (0,1) values from the Acceptance-Rejection method
// @ MTH9821 Homework1 Group6 

#include "ARM.h"
#include <random>

#include "LCG.h"

#include <cmath>

double ARM::standard_normal() {


    double u1, u2, u3, x;

    // looping while the can accept model
    while (true) {

        // generate three numbers from uniform distribution on [0,1]
        u1 = LCG::gen();
        u2 = LCG::gen();
        u3 = LCG::gen();

        // generate number from an exponential with param 1
        x = -std::log(u1);

        // checking if we can accept this number
        if (u2 < std::exp(-(x - 1.) * (x - 1.) / 2.)) {
            // accept
            break;
        };
    };
    // swhitch sign with .5 probobility
    if (u3 < .5) x = -x;
    return x;
};


// generate rv and count the number of uniform variable used
std::tuple<int, double> ARM::standard_normal_count() {
    int num = 0;
    // looping while the can accept model
    double u1, u2, u3, x;

    // looping while the can accept model
    while (true) {

        // generate three numbers from uniform distribution on [0,1]
        u1 = LCG::gen();
        u2 = LCG::gen();
        u3 = LCG::gen();
        num += 3;
        // generate number from an exponential with param 1
        x = -std::log(u1);

        // checking if we can accept this number
        if (u2 < std::exp(-(x - 1.) * (x - 1.) / 2.)) {
            // accept
            break;
        }
    }
    // swhitch sign with .5 probobility
    if (u3 < .5) 
        x = -x;
    return std::make_tuple(num, x);


    //std::random_device rd{};
    //std::mt19937 gen{ rd() };

    //// values near the mean are the most likely
    //// standard deviation affects the dispersion of generated values from the mean
    //std::normal_distribution<> d{ 0,1 };
    //double x = d(gen);
    //return std::make_tuple(3, x);
}