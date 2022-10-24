//
//  FiniteDifferencePricer.cpp
//  MTH_9821_HW7_Group6
//
//  Created by 王明森 on 2022/10/21.
//

#include "FiniteDifferencePricer.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

#include <iostream>

double FiniteDifferencePricer::alpha_temp_ = 0.45;
std::size_t FiniteDifferencePricer::M_init_ = 4;

FiniteDifferencePricer::FiniteDifferencePricer(double S0, double K, double T, double sigma, double r, double q) : S0_(S0), K_(K), T_(T), sigma_(sigma), q_(q), r_(r) {
    
    //  1. Computational domain (part 1)
    double sigma2 = sigma * sigma;
    
    tau_final_ = T * sigma2 / 2.;
    
    x_l_ = std::log(S0 / K) + (r - q - sigma2 / 2.) * T - 3. * sigma * std::sqrt(T);
    x_r_ = std::log(S0 / K) + (r - q - sigma2 / 2.) * T + 3. * sigma * std::sqrt(T);
    
    
    
    // Heat equation transformation coefficients
    a_ = (r - q) / sigma2 - .5;
    b_ = (a_ + 1.) * (a_ + 1.) + 2. * q / sigma2;
    
}

void FiniteDifferencePricer::EuroPut_advance(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh) const {
    
    std::vector<double> new_u_mesh;
    
    // Left boundary
    new_u_mesh.push_back(boundary_x_l_(tau));
    
    // Middle values
    for (std::size_t pos = 1; pos < x_mesh.size() - 1; pos++) {
        new_u_mesh.push_back(alpha * u_mesh[pos - 1] + (1. - 2. * alpha) * u_mesh[pos] + alpha * u_mesh[pos + 1]);
    }
    
    // Right boundary
    new_u_mesh.push_back(boundary_x_r_(tau));
    
    u_mesh = std::move(new_u_mesh);
}

void FiniteDifferencePricer::PrintVector(const std::vector<double>& vec) const {
    for (auto elem : vec) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}

std::vector<double> FiniteDifferencePricer::EuroPut(std::size_t M) {
    
    // 2. Boundary conditions
    boundary_tau_0_ = [=](double x)->double {
        // K * exp(ax) * max(1 - exp(x), 0)
        if (x < 0.) {
            return K_ * std::exp(a_ * x) * (1. - std::exp(x));
        } else {
            return 0.;
        }
    };
    
    boundary_x_l_ = [=](double tau)->double {
        // K * exp(ax_{left} + b\tau)(exp(-2r\tau / sigma^2) - exp(x_{left} - 2q\tau / sigma^2))
        return K_ * std::exp(a_ * x_l_ + b_ * tau) * (std::exp(-2. * r_ * tau / (sigma_ * sigma_)) - std::exp(x_l_ - 2. * q_ * tau / (sigma_ * sigma_)));
    };
    
    boundary_x_r_ = [](double tau)->double {
        return 0.;
    };
    
    return this->EuroPut_impl(M);
    
}

std::vector<double> FiniteDifferencePricer::EuroPut_impl(std::size_t M) const {
    
    std::vector<double> res;
    
    // 1. Computational domain (part 2)
    // M: number of time intervals on the tau-axis
    // N: number of x intervals on the x-axis
    std::size_t N;
    // Interval lengths
    double dtau, dx, alpha;
    std::tie(N, dtau, dx, alpha) = this->DomainParams(M);
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    std::size_t interval_i;
    std::tie(x_mesh, interval_i) = this->BuildMesh(N, dx);
    
    std::vector<std::vector<double>> u_meshes = this->FiniteDifference(alpha, x_mesh, M, dtau);
    std::vector<double> u_mesh = u_meshes[0];
    std::vector<double> u_mesh_prev = u_meshes[1];
    
    // 4. Pointwise convergence
    std::vector<double> approximations = this->Approximate(x_mesh, u_mesh, dx);
    res.insert(res.end(), approximations.cbegin(), approximations.cend());
    double V_approx = approximations.front();
    
    // 5. RMS error
    double error_RMS = this->RMS(x_mesh, u_mesh);
    res.push_back(error_RMS);
    
    // 6. Greeks
    std::vector<double> Greeks = this->Greeks(interval_i, x_mesh, u_mesh, u_mesh_prev, dtau, V_approx);
    res.insert(res.end(), Greeks.cbegin(), Greeks.cend());
    
    return res;
}

std::tuple<std::size_t, double, double, double> FiniteDifferencePricer::DomainParams(std::size_t M) const {
    
    // More finite difference parameters
    double dtau = tau_final_ / M;
    // N: number of x intervals on the x-axis
    std::size_t N = std::floor((x_r_ - x_l_) / std::sqrt(dtau / alpha_temp_));
    double dx = (x_r_ - x_l_) / N;
    double alpha = dtau / (dx * dx);
    
    return std::make_tuple(N, dtau, dx, alpha);
    
}

std::pair<std::vector<double>, std::size_t> FiniteDifferencePricer::BuildMesh(std::size_t N, double dx) const {
    
    // Fill x mesh
    std::vector<double> x_mesh({x_l_});
    for (std::size_t i = 1; i < N; i++) {
        x_mesh.push_back(x_mesh.back() + dx);
    }
    
    // Find actual x on x_mesh
    double x_compute = std::log(S0_ / K_);
    auto x_large_it = std::upper_bound(x_mesh.cbegin(), x_mesh.cend(), x_compute);
    std::size_t interval_i = std::distance(x_mesh.cbegin(), x_large_it) - 1;
    
    return std::make_pair(x_mesh, interval_i);
}

std::vector<std::vector<double>> FiniteDifferencePricer::FiniteDifference(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const {
    
    // Fill u mesh with boundary condition
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), boundary_tau_0_);
    //this->PrintVector(u_mesh);
    // Advance M-1 times
    for (std::size_t i = 1; i < M; i++) {
        double curr_tau = dtau * i;
        this->EuroPut_advance(curr_tau, alpha, x_mesh, u_mesh);
        //this->PrintVector(u_mesh);
    }
    
    // Record the (M-1)th u_mesh for theta calculation
    std::vector<double> u_mesh_prev = u_mesh;
    this->EuroPut_advance(tau_final_, alpha, x_mesh, u_mesh);
    //this->PrintVector(u_mesh);
    
    return std::vector<std::vector<double>>({u_mesh, u_mesh_prev});
}

std::vector<double> FiniteDifferencePricer::Approximate(const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, double dx) const {
    
    std::vector<double> approximations;
    
    // Find the interval containing x_compute
    double x_compute = std::log(S0_ / K_);
    auto x_large_it = std::upper_bound(x_mesh.cbegin(), x_mesh.cend(), x_compute);
    std::size_t interval_i = std::distance(x_mesh.cbegin(), x_large_it) - 1;
    
    // Approximation method 1
    double S_small = K_ * std::exp(x_mesh[interval_i]);
    double S_large = K_ * std::exp(x_mesh[interval_i + 1]);
    double V_small = std::exp(-a_ * (x_mesh[interval_i]) - b_ * tau_final_) * u_mesh[interval_i];
    double V_large = std::exp(-a_ * (x_mesh[interval_i + 1]) - b_ * tau_final_) * u_mesh[interval_i + 1];
    double V_approx = ((S_large - S0_) * V_small + (S0_ - S_small) * V_large) / (S_large - S_small);
    approximations.push_back(V_approx);
    
    // Approximation method 2
    double u_approx = ((x_mesh[interval_i + 1] - x_compute) * u_mesh[interval_i] + (x_compute - x_mesh[interval_i]) * u_mesh[interval_i + 1]) / (dx);
    double V_approx_2 = std::exp(-a_ * x_compute - b_ * tau_final_) * u_approx;
    approximations.push_back(V_approx_2);
    
    return approximations;
}

double FiniteDifferencePricer::RMS(const std::vector<double>& x_mesh, const std::vector<double>& u_mesh) const {
    
    // Calculate vector of FD and BS option values
    auto V_mesh_approx_gen = [&](double x, double u)->double {
        return std::exp(-a_ * x - b_ * tau_final_) * u;
    };
    auto V_mesh_exact_gen = [&](double x)->double {
        EuropeanOption option(0., K_ * std::exp(x), K_, T_, sigma_, r_, q_);
        return option.Put();
    };
    
    std::vector<double> V_mesh_approx(x_mesh.size());
    std::vector<double> V_mesh_exact(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.cbegin(), V_mesh_approx.begin(), V_mesh_approx_gen);
    std::transform(x_mesh.cbegin(), x_mesh.cend(), V_mesh_exact.begin(), V_mesh_exact_gen);
    
    // Find RMS
    double error_sq = 0.;
    int error_count = 0;
    auto V_mesh_approx_it = V_mesh_approx.cbegin();
    auto V_mesh_exact_it = V_mesh_exact.cbegin();
    while (V_mesh_exact_it != V_mesh_exact.cend()) {
        double V_BS = *V_mesh_exact_it;
        double V_FD = *V_mesh_approx_it;
        if (V_BS > 0.00001 * S0_) {
            error_count++;
            error_sq += (V_BS - V_FD) * (V_BS - V_FD) / (V_BS * V_BS);
        }
        V_mesh_approx_it++;
        V_mesh_exact_it++;
    }
    double error_RMS = std::sqrt(error_sq / error_count);
    
    return error_RMS;
}

std::vector<double> FiniteDifferencePricer::Greeks(std::size_t interval_i, const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, std::vector<double>& u_mesh_prev, double dtau, double V_approx) const {
    
    // DELTA & GAMMA
    // Get S and V of interest
    double S_smaller = K_ * std::exp(x_mesh[interval_i - 1]);
    double S_small = K_ * std::exp(x_mesh[interval_i]);
    double S_large = K_ * std::exp(x_mesh[interval_i + 1]);
    double S_larger = K_ * std::exp(x_mesh[interval_i + 2]);
    
    double V_smaller = std::exp(-a_ * (x_mesh[interval_i - 1]) - b_ * tau_final_) * u_mesh[interval_i - 1];
    double V_small = std::exp(-a_ * (x_mesh[interval_i]) - b_ * tau_final_) * u_mesh[interval_i];
    double V_large = std::exp(-a_ * (x_mesh[interval_i + 1]) - b_ * tau_final_) * u_mesh[interval_i + 1];
    double V_larger = std::exp(-a_ * (x_mesh[interval_i + 2]) - b_ * tau_final_) * u_mesh[interval_i + 2];
    
    // Find delta and gamma
    double delta = (V_large - V_small) / (S_large - S_small);
    double gamma = ((V_larger - V_large) / (S_larger - S_large) - (V_small - V_smaller) / (S_small - S_smaller)) / (((S_larger + S_large) / 2.) - ((S_small + S_smaller) / 2.));
    
    // THETA
    // Get dt from dtau
    double dt = 2. * dtau / (sigma_ * sigma_);
    
    // Get V at t = dt
    double V_small_prev = std::exp(-a_ * (x_mesh[interval_i]) - b_ * (tau_final_ - dtau)) * u_mesh_prev[interval_i];
    double V_large_prev = std::exp(-a_ * (x_mesh[interval_i + 1]) - b_ * (tau_final_ - dtau)) * u_mesh_prev[interval_i + 1];
    double V_approx_prev = ((S_large - S0_) * V_small_prev + (S0_ - S_small) * V_large_prev) / (S_large - S_small);
    double theta = (V_approx_prev - V_approx) / dt;
    
    return std::vector<double>({delta, gamma, theta});
}


std::vector<double> FiniteDifferencePricer::AmericanPut(std::size_t M) {
    
    // 2. Boundary conditions
    boundary_tau_0_ = [=](double x)->double {
        // K * exp(ax) * max(1 - exp(x), 0)
        if (x < 0.) {
            return K_ * std::exp(a_ * x) * (1. - std::exp(x));
        } else {
            return 0.;
        }
    };
    
    boundary_x_l_ = [=](double tau)->double {
        // K * exp(ax_{left} + b\tau)(exp(-2r\tau / sigma^2) - exp(x_{left} - 2q\tau / sigma^2))
        return K_ * std::exp(a_ * x_l_ + b_ * tau) * (std::exp(-2. * r_ * tau / (sigma_ * sigma_)) - std::exp(x_l_ - 2. * q_ * tau / (sigma_ * sigma_)));
    };
    
    boundary_x_r_ = [](double tau)->double {
        return 0.;
    };
    
    return this->AmericanPut_impl(M);
}

std::vector<double> FiniteDifferencePricer::AmericanPut_impl(std::size_t M) const {
    
    std::vector<double> res;
    
    // 1. Computational domain (part 2)
    // M: number of time intervals on the tau-axis
    // N: number of x intervals on the x-axis
    std::size_t N;
    // Interval lengths
    double dtau, dx, alpha;
    std::tie(N, dtau, dx, alpha) = this->DomainParams(M);
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    std::size_t interval_i;
    std::tie(x_mesh, interval_i) = this->BuildMesh(N, dx);
    
    std::vector<std::vector<double>> u_meshes = this->FiniteDifference_AmeriPut(alpha, x_mesh, M, dtau);
    std::vector<double> u_mesh = u_meshes[0];
    std::vector<double> u_mesh_prev = u_meshes[1];
    
    // 4. Pointwise convergence
    std::vector<double> approximations = this->Approximate(x_mesh, u_mesh, dx);
    res.insert(res.end(), approximations.cbegin(), approximations.cend());
    double V_approx = approximations.front();
    
    // 6. Greeks
    std::vector<double> Greeks = this->Greeks(interval_i, x_mesh, u_mesh, u_mesh_prev, dtau, V_approx);
    res.insert(res.end(), Greeks.cbegin(), Greeks.cend());
    
    // 7. Variance reduction
    res.push_back(this->VarianceReduction_AmeriPut(V_approx, M));
    
    
    return res;
    
}

std::vector<std::vector<double>> FiniteDifferencePricer::FiniteDifference_AmeriPut(double alpha, const std::vector<double>& x_mesh, std::size_t M, double dtau) const {
    
    // Fill u mesh with boundary condition
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), boundary_tau_0_);
    //this->PrintVector(u_mesh);
    
    // Advance M-1 times
    for (std::size_t i = 1; i < M; i++) {
        double curr_tau = dtau * i;
        this->AmeriPut_advance(curr_tau, alpha, x_mesh, u_mesh);
        //this->PrintVector(u_mesh);
    }
    
    // Record the (M-1)th u_mesh for theta calculation
    std::vector<double> u_mesh_prev = u_mesh;
    this->AmeriPut_advance(tau_final_, alpha, x_mesh, u_mesh);
    //this->PrintVector(u_mesh);

    return std::vector<std::vector<double>>({u_mesh, u_mesh_prev});
}

void FiniteDifferencePricer::AmeriPut_advance(double tau, double alpha, const std::vector<double>& x_mesh, std::vector<double>& u_mesh) const {
    
    std::vector<double> new_u_mesh;
    
    // Left boundary
    new_u_mesh.push_back(boundary_x_l_(tau));
    
    // Middle values
    for (std::size_t pos = 1; pos < x_mesh.size() - 1; pos++) {
        // Get the corresponding European option's value
        double euro_val = alpha * u_mesh[pos - 1] + (1. - 2. * alpha) * u_mesh[pos] + alpha * u_mesh[pos + 1];
        
        // Find early exercise
        double early_ex_premium = 0.;
        
        if (x_mesh[pos] < 0.) {
            early_ex_premium = K_ * std::exp(a_ * x_mesh[pos] + b_ * tau) * (1. - std::exp(x_mesh[pos]));
        }
        
        // Compare and add to mesh
        new_u_mesh.push_back(std::max(euro_val, early_ex_premium));
    }
    
    // Right boundary
    new_u_mesh.push_back(boundary_x_r_(tau));
    
    u_mesh = std::move(new_u_mesh);
}


double FiniteDifferencePricer::VarianceReduction_AmeriPut(double V_approx, std::size_t M) const {
    
    // Get corresponding European option value for BS and FD
    FiniteDifferencePricer FDPricer(S0_, K_, T_, sigma_, r_, q_);
    EuropeanOption BSPricer(0., S0_, K_, T_, sigma_, r_, q_);
    
    double FDEuro = FDPricer.EuroPut(M).front();
    double BSEuro = BSPricer.Put();
    
    // Adjust approximation by the pointwise difference
    return V_approx + BSEuro - FDEuro;
    
}


std::array<std::vector<double>, 2> FiniteDifferencePricer::AmericanPut_EarlyExDomain(std::size_t M) {
    
    // 2. Boundary conditions
    boundary_tau_0_ = [=](double x)->double {
        // K * exp(ax) * max(1 - exp(x), 0)
        if (x < 0.) {
            return K_ * std::exp(a_ * x) * (1. - std::exp(x));
        } else {
            return 0.;
        }
    };
    
    boundary_x_l_ = [=](double tau)->double {
        // K * exp(ax_{left} + b\tau)(exp(-2r\tau / sigma^2) - exp(x_{left} - 2q\tau / sigma^2))
        return K_ * std::exp(a_ * x_l_ + b_ * tau) * (std::exp(-2. * r_ * tau / (sigma_ * sigma_)) - std::exp(x_l_ - 2. * q_ * tau / (sigma_ * sigma_)));
    };
    
    boundary_x_r_ = [](double tau)->double {
        return 0.;
    };
    
    // 1. Computational domain (part 2)
    // M: number of time intervals on the tau-axis
    // N: number of x intervals on the x-axis
    std::size_t N;
    // Interval lengths
    double dtau, dx, alpha;
    std::tie(N, dtau, dx, alpha) = this->DomainParams(M);
    
    // 3. Finite difference scheme
    // Fill x mesh
    std::vector<double> x_mesh;
    std::size_t interval_i;
    std::tie(x_mesh, interval_i) = this->BuildMesh(N, dx);
    
    // MODIFIED FINITE DIFFERENCE
    // Fill u mesh with boundary condition
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), boundary_tau_0_);
    
    // Advance M times while recording early exercise position
    std::vector<double> Sopt;
    for (std::size_t i = 1; i <= M; i++) {
        // Record maximum exercise position
        std::size_t Nopt = 0;
        
        double tau = dtau * i;
        
        std::vector<double> new_u_mesh;
        
        // Left boundary
        new_u_mesh.push_back(boundary_x_l_(tau));
        
        // Middle values
        for (std::size_t pos = 1; pos < x_mesh.size() - 1; pos++) {
            // Get the corresponding European option's value
            double euro_val = alpha * u_mesh[pos - 1] + (1. - 2. * alpha) * u_mesh[pos] + alpha * u_mesh[pos + 1];
            
            // Find early exercise
            double early_ex_premium = 0.;
            
            if (x_mesh[pos] < 0.) {
                early_ex_premium = K_ * std::exp(a_ * x_mesh[pos] + b_ * tau) * (1. - std::exp(x_mesh[pos]));
            }
            
            // Compare and add to mesh
            if (early_ex_premium >= euro_val) {
                new_u_mesh.push_back(early_ex_premium);
                if (early_ex_premium > 0) Nopt = pos;
            } else {
                new_u_mesh.push_back(euro_val);
            }
        }
        double S_small = K_ * std::exp(x_mesh[Nopt]);
        double S_large = K_ * std::exp(x_mesh[Nopt + 1]);
        Sopt.push_back((S_small + S_large) / 2.);
        
        // Right boundary
        new_u_mesh.push_back(boundary_x_r_(tau));
        
        u_mesh = std::move(new_u_mesh);
    }
    
    // Get array of t
    std::vector<double> t_mesh(16);
    std::iota(t_mesh.begin(), t_mesh.end(), 1.);
    auto convert_to_t = [&](double m)->double {
        return T_ - (2. * m * dtau) / (sigma_ * sigma_);
    };
    std::transform(t_mesh.begin(), t_mesh.end(), t_mesh.begin(), convert_to_t);
    
//    this->PrintVector(x_mesh);
    return std::array<std::vector<double>, 2>({t_mesh, Sopt});
}
