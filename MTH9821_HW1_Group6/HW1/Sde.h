// ISde: class that evaluates Stochastic differential equations
// taken out of the Advanced C++ course
// @ MTH9821 Homework1 Group6 

#ifndef Sde_h
#define Sde_h

class ISde {
    // Abstract class for all standard one-factor SDEs
    // dX = mu(X,t)dt + sigma(X,t)dW, given X(0)
    // Convection-diffusion
    
    double S0_;     // Initial condition (Asset price at t = 0)
    double expiry_; // Expiry
    
public:
    ISde(double S0, double expiry);
    virtual ~ISde() = default;
    
    virtual double Convection(double x, double t) const = 0;    // mu
    virtual double Diffusion(double x, double t) const = 0;     // sigma
    
    // Some extra functions
    virtual double ConvectionCorrected(double x, double t, double B) const = 0;
    virtual double DiffusionDerivative(double x, double t) const = 0;
    
    // Get / set initial condition
    virtual double S0() const final;
    virtual void S0(double x0) final;
    
    // Get / set expiry time
    virtual double Expiry() const final;  // T
    virtual void Expiry(double expiry) final;
};

class GBM : public ISde {
    // Geometric Brownian motion
private:
    double r_;      // Riskless rate
    double vol_;    // Volatility (constant)
    double d_;      // Dividend
    
    
public:
    GBM(double riskless_rate, double volatility, double dividend_rate, double initial_condition, double expiry);
    virtual ~GBM() override final = default;
    
    // Get coefficients
    virtual double Convection(double x, double t) const override final;
    virtual double Diffusion(double x, double t) const override final;
    virtual double ConvectionCorrected(double x, double t, double B) const override final;
    virtual double DiffusionDerivative(double x, double t) const override final;
};

class CEV : public ISde {
    // Constant elasticity of variance
private:
    double r_;      // Riskless rate
    double vol_;    // Volatility (constant)
    double d_;      // Dividend
    double beta_;   // CEV parameter
    
public:
    CEV(double riskless_rate, double volatility, double dividend_rate, double initial_condition, double expiry, double beta);
    virtual ~CEV() override final = default;
    
    // Get coefficients
    virtual double Convection(double x, double t) const override final;
    virtual double Diffusion(double x, double t) const override final;
    virtual double ConvectionCorrected(double x, double t, double B) const override final;
    virtual double DiffusionDerivative(double x, double t) const override final;
};

#endif /* Sde_h */
