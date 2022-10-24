//
//  FiniteDifferencePricer.hpp
//  MTH_9821_HW7_Group6
//
//  Created by 王明森 on 2022/10/21.
//

#ifndef FiniteDifferencePricer_hpp
#define FiniteDifferencePricer_hpp

#include "EuropeanOption.hpp"
#include <functional>
#include <vector>
#include <tuple>
#include <array>

class FiniteDifferencePricer {
private:
    // Option data
    double S0_;
    double K_;
    double T_;
    double sigma_;
    double r_;
    double q_;
    
    // Heat equation transformation coefficients
    double a_;
    double b_;
    
    // Finite difference hyperparameters
    static double alpha_temp_;
    static std::size_t M_init_;
    
    // Finite difference parameters
    double tau_final_;
    double x_l_;
    double x_r_;
//    double dt_;
//    double dx_;
//    std::size_t N_;
//    double alpha_;
    
    // Boundary conditions
    std::function<double (double)> boundary_tau_0_;
    std::function<double (double)> boundary_x_l_;
    std::function<double (double)> boundary_x_r_;
    
    
    // 1. Computational domain
    std::tuple<std::size_t, double, double, double> DomainParams(std::size_t M) const;
    
    std::pair<std::vector<double>, std::size_t> BuildMesh(std::size_t N, double dx) const;
    
    // 3. Finite difference scheme
    std::vector<std::vector<double>> FiniteDifference(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const;
    
    std::vector<std::vector<double>> FiniteDifference_AmeriPut(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const;
    
    // Advance time (modify u-mesh in place)
    void EuroPut_advance(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh) const;
    void AmeriPut_advance(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh) const;
    
    // 4. Pointwise convergence
    std::vector<double> Approximate(const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, double dx) const;
    
    // 5. RMS error
    double RMS(const std::vector<double>& x_mesh, const std::vector<double>& u_mesh) const;
    
    // 6. Greeks
    std::vector<double> Greeks(std::size_t interval_i, const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, std::vector<double>& u_mesh_prev, double dtau, double V_approx) const;
    
    // 7. Variance reduction
    double VarianceReduction_AmeriPut(double V_approx, std::size_t M) const;
    
public:
    FiniteDifferencePricer(double S0, double K, double T, double sigma, double r, double q);
    ~FiniteDifferencePricer() = default;
    void PrintVector(const std::vector<double>& vec) const;
    
    std::vector<double> EuroPut(std::size_t M);
    std::vector<double> EuroPut_impl(std::size_t M) const;
    
    std::vector<double> AmericanPut(std::size_t M);
    std::vector<double> AmericanPut_impl(std::size_t M) const;
    
    std::array<std::vector<double>, 2> AmericanPut_EarlyExDomain(std::size_t M);
};

#endif /* FiniteDifferencePricer_hpp */
