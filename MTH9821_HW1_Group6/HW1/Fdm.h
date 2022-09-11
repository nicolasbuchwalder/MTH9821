// FDM: multiple classes for approximating SDES
// // taken out of the Advanced C++ course
// @ MTH9821 Homework1 Group6 

#ifndef Fdm_h
#define Fdm_h

#include "Sde.h"
#include <memory>

class IFdm {
    // Abstract class for one-step FDM to approximate one-factor SDEs
public:
    virtual ~IFdm() = default;
    
    // Get / set underlying SDE
    virtual std::shared_ptr<ISde> SDE() const = 0;
    virtual void SDE(const std::shared_ptr<ISde>& sde) = 0;
    
    // Advance solution from t to t+dt
    virtual double advance(double x, unsigned step, double wiener_increment) const = 0;
};

class FdmBase : public IFdm {
    // Base class of FDM (still abstract)
protected:
    std::shared_ptr<ISde> sde_; // Underlying SDE
    unsigned NT_;               // # of time steps
    double dt_;                 // Time increment
    double dt_sqrt_;            // Square root of time increment
    
    
public:
    FdmBase(const std::shared_ptr<ISde>& sde, unsigned NT);
    virtual ~FdmBase() override = default;
    
    // Get / set underlying SDE
    virtual std::shared_ptr<ISde> SDE() const override final;
    virtual void SDE(const std::shared_ptr<ISde>& sde) override final;
};

class Euler : public FdmBase {
    // One-factor explicit Euler
public:
    Euler(const std::shared_ptr<ISde>& sde, unsigned NT);
    virtual ~Euler() override final = default;
    
    virtual double advance(double x, unsigned step, double wiener_increment) const override final;
};

class Multiplicative : public FdmBase {
    // Multiplicative method
private:
    double sig_;    // Diffusion (volatility)
    double mu_;     // Convection
    
public:
    Multiplicative(const std::shared_ptr<ISde>& sde, unsigned NT, double sig, double mu);
    virtual ~Multiplicative() override final = default;
    
    virtual double advance(double x, unsigned step, double wiener_increment) const override final;
};

class Exact : public FdmBase {
    // Compute exact value at t + dt
private:
    double S0_;     // Initial condition (initial value of the underlying asset)
    double sig_;    // Diffusion (volatility)
    double mu_;     // Convection
    
public:
    Exact(const std::shared_ptr<ISde>& sde, unsigned NT, double sig, double mu);
    virtual ~Exact() override final = default;
    
    virtual double advance(double x, unsigned step, double wiener_increment) const override final;
};

class Milstein : public FdmBase {
    // Milstein method
    // https://en.wikipedia.org/wiki/Milstein_method#Intuitive_derivation
public:
    Milstein(const std::shared_ptr<ISde>& sde, unsigned NT);
    virtual ~Milstein() override final = default;
    
    virtual double advance(double x, unsigned step, double wiener_increment) const override final;
};

class PredictorCorrector : public FdmBase {
    // Predictor-Corrector method
    // https://en.wikipedia.org/wiki/Predictor–corrector_method
    
private:
    double a_;
    double b_;
    
public:
    PredictorCorrector(const std::shared_ptr<ISde>& sde, unsigned NT, double a, double b);
    virtual ~PredictorCorrector() override final = default;
    
    virtual double advance(double x, unsigned step, double wiener_increment) const override final;
};

class ModifiedPredictorCorrector : public FdmBase {
    // Modified predictor-Corrector method
    // https://en.wikipedia.org/wiki/Predictor–corrector_method
private:
    double a_;
    double b_;
    
public:
    ModifiedPredictorCorrector(const std::shared_ptr<ISde>& sde, unsigned NT, double a, double b);
    virtual ~ModifiedPredictorCorrector() override final = default;
    
    virtual double advance(double x, unsigned step, double wiener_increment) const override final;
};


#endif /* Fdm_h */
