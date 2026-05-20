#include "CoulCorrCalc.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>

// Constructor
CoulCorrCalc::CoulCorrCalc()
{
  HypCalculatorInstance = new HypCalculator();
  NMaxIter = 3;
  epsTolerance = 1e-2;
  Mc2 = Mass_Pi*1000.;
  betaT = 0;
}

// Destructor
CoulCorrCalc::~CoulCorrCalc()
{
  delete HypCalculatorInstance;
}

// Set integration properties
void CoulCorrCalc::SetIntegrationProperties(int _NMaxIter, double _epsTolerance)
{
  NMaxIter = _NMaxIter;
  epsTolerance = _epsTolerance;
  
}

// Integrand of the spherical averaging of f_s, with y=cos(theta), with Rcc in LCMS
double CoulCorrCalc::f_s_y_integrand(const double y, const double q, const double RccLCMS, const double alpha)
{
  return exp(-0.5 * pow(fabs(q*RccLCMS / (HBARC*1000.)), alpha) * pow(1.+betaT*betaT/(1.-betaT*betaT)*y*y, alpha/2.));
}

// No-Coulomb auxiliary correlation function term, with the pair-radius Rcc
double CoulCorrCalc::f_s(const double q, const double Rcc, const double alpha)
{
  if(betaT==0)
    return exp(-0.5 * pow(fabs(q*Rcc / (HBARC*1000.)), alpha)); // Rcc is in PCMS, spherical symmetry assumed in PCMS
  else
  {
    double error = 0.;
    double lowerlim = -1.;
    double upperlim = 1.;
    //auto func = bind(&CoulCorrCalc::f_s_y_integrand, this, placeholders::_1, q, Rcc, alpha); // Rcc is in LCMS, spherical symmetry assumed in LCMS
    auto func = [this, q, Rcc, alpha](double y) { return f_s_y_integrand(y, q, Rcc, alpha); }; // Rcc is in LCMS, spherical symmetry assumed in LCMS
    double result = 0.5*boost::math::quadrature::gauss_kronrod<double, NGaussKronrod>::integrate(func, lowerlim, upperlim, NMaxIter, epsTolerance, &error);
    return result;
  }
}

// Exact result for A_1,s
double CoulCorrCalc::A_1_s_wo_int(const double x, const double k, const double Rcc, const double alpha, const double eta)
{
  NFuncCalls++; //This is just for monitoring, can be commented out
  double fs1_1 = (f_s(2.*k*x, Rcc, alpha) - f_s(0., Rcc, alpha)) / x;
  double fs1_2 = (f_s(2.*k/x, Rcc, alpha) - f_s(0., Rcc, alpha)) / x;
  complex<double> func_1 = pow(1. + 1./x, 2.*eta*I) * HypCalculatorInstance->Fplus(1. / (x*x));
  complex<double> func_2 = pow(1. + x, 2.*eta*I) *  HypCalculatorInstance->Fplus(x*x);
  return -2./eta * imag(fs1_1 * func_1 + fs1_2 * func_2);
}

// Exact result for A_2,s
double CoulCorrCalc::A_2_s_wo_int(const double x, const double k, const double Rcc, const double alpha, const double eta)
{
  NFuncCalls++; //This is just for monitoring, can be commented out
  double func_1 = sin(eta * log((1. + x)/(1. - x))) / (x * (x + 1.));
  double func_2 = x * x * exp(M_PI * eta) * (f_s(2. * k * x, Rcc, alpha) - f_s(2. * k, Rcc, alpha)) / (1. - x);
  double func_3 = (f_s(2. * k / x, Rcc, alpha) - f_s(2. * k, Rcc, alpha)) / (1. - x);
  return (2. / eta ) * (func_1 * func_2 - func_1 * func_3);
}

// The A1 numerical integral, with boost::math::quadrature::gauss_kronrod
double CoulCorrCalc::A1_int(const double k, const double Rcc, const double alpha, const double eta)
{
  HypCalculatorInstance->initialize_eta(eta);

  double error = 0;
  double lowerlim = 0;
  double upperlim = 1;
  //auto func = bind(&CoulCorrCalc::A_1_s_wo_int, this, placeholders::_1, k, Rcc, alpha, eta);
  auto func = [this, k, Rcc, alpha, eta](double x) { return A_1_s_wo_int(x, k, Rcc, alpha, eta); };

  double result = boost::math::quadrature::gauss_kronrod<double, NGaussKronrod>::integrate(func, lowerlim, upperlim, NMaxIter, epsTolerance, &error);
  return result;
}

// The A2 numerical integral, with boost::math::quadrature::gauss_kronrod
double CoulCorrCalc::A2_int(const double k, const double Rcc, const double alpha, const double eta)
{
  HypCalculatorInstance->initialize_eta(eta);

  double error = 0;
  double lowerlim = 0;
  double upperlim = 1;
  //auto func = bind(&CoulCorrCalc::A_2_s_wo_int, this, placeholders::_1, k, Rcc, alpha, eta);
  auto func = [this, k, Rcc, alpha, eta](double x) { return A_2_s_wo_int(x, k, Rcc, alpha, eta); };

  double result = boost::math::quadrature::gauss_kronrod<double, NGaussKronrod>::integrate(func, lowerlim, upperlim, NMaxIter, epsTolerance, &error);
  return result;
}

// Set the particle mass, argument in GeV/c^2
void CoulCorrCalc::SetParticleMass(const double mc2)
{
  Mc2 = mc2*1000; // Convert to MeV
}

// Set betaT
void CoulCorrCalc::SetBetaT(const double betat)
{
  betaT = betat;
}

// Calculating the Sommerfeld parameter eta
const double CoulCorrCalc::calc_eta(const double k)
{
  return FINESTRUCTURE_CONSTANT * Mc2 * 0.5 / k;
}

// Correlation function, for lambda=1, with Coulomb
double CoulCorrCalc::FullCorrFuncValue(const double alpha, const double R, const double Q)
{
  NFuncCalls = 0;
  double k = Q*500; //k = Q/2, but in MeV here
  double Rcc = R*pow(2.,1./alpha);
  double eta = calc_eta(k);
  double Gamow = (2. * M_PI * eta) / (exp(2. * M_PI * eta) - 1.);
  double A1 = A1_int(k, Rcc, alpha, eta);
  double A2 = A2_int(k, Rcc, alpha, eta);
  double value = Gamow * (1. + f_s(2.*k, Rcc, alpha) + (eta / M_PI) * (A1 + A2));
  return value;
}

// Correlation function, incorporating the lambda value, via the Bowler-Sinyikov formula
double CoulCorrCalc::FullCorrFuncValueLambda(const double alpha, const double R, double lambda, const double Q)
{
  return 1 - lambda + lambda * FullCorrFuncValue(alpha, R, Q);
}

// Correlation function, without Coulomb, with lambda
double CoulCorrCalc::PureCorrFuncValueLambda(const double alpha, const double R, double lambda, const double Q)
{
  double Rcc = R*pow(2.,1./alpha);
  return 1.0 + lambda * f_s(Q*1000, Rcc, alpha);
}

// Coulomb correction, without the lambda value
double CoulCorrCalc::CoulCorrValue(const double alpha, const double R, const double Q)
{
  double C0value = PureCorrFuncValueLambda(alpha, R, 1.0, Q);
  double CCvalue = FullCorrFuncValue(alpha, R, Q);
  return CCvalue/C0value;
}

