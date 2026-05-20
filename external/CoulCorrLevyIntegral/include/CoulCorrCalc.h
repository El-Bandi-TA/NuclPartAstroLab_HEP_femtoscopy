#ifndef _CoulCorrCalc_h_
#define _CoulCorrCalc_h_

//#include <algorithm>
//#include <functional>
//#include <iostream>
//#include <cmath>
//#include <limits>
#include "basics.h"
#include "HypCalculator.h"
#include "CoulCorrCalc.h"

class CoulCorrCalc
{
 public:
  CoulCorrCalc();
  ~CoulCorrCalc();
  
  // Correlation function, for lambda=1, with Coulomb, Q in GeV/c, R in fm
  double FullCorrFuncValue(const double alpha, const double R, const double Q);
  
  // Correlation function, incorporating the lambda value, via the Bowler-Sinyikov formula
  double FullCorrFuncValueLambda(const double alpha, const double R, double lambda, const double Q);
  
  // Correlation function, without Coulomb, with lambda
  double PureCorrFuncValueLambda(const double alpha, const double R, double lambda, const double Q);
  
  // Coulomb correction, without the lambda value
  double CoulCorrValue(const double alpha, const double R, const double Q);
  
  // Set integration properties
  void SetIntegrationProperties(int NMaxIter, double epsTolerance);
  
  // Get number of function calls in last calculation
  int GetNFuncCalls() { return NFuncCalls; }

  // Set the particle mass (the actual mass, not the reduced one), in GeV/c^2
  void SetParticleMass(const double mc2);

  // Set betaT = KT/sqrt(mT^2+KT^2)
  void SetBetaT(const double betat);

private:
  // Private functions, explained in the source code
  inline double f_s_y_integrand(const double y, const double q, const double RLCMS, const double alpha);
  double f_s(const double q, const double R, const double alpha);
  inline double A_1_s_wo_int(const double x, const double k, const double R, const double alpha, const double eta);
  inline double A_2_s_wo_int(const double x, const double k, const double R, const double alpha, const double eta);
  double A1_int(const double k, const double R, const double alpha, const double eta);
  double A2_int(const double k, const double R, const double alpha, const double eta);
  const double calc_eta(const double k); // k: MeV/c

  // Physical properties
  double Mc2; // Mass of the particle in MeV/c^2
  double betaT; // betaT between PCMS and LCMS

  // An instance of the hypergeometric 2F1 calculator
  HypCalculator* HypCalculatorInstance;
  
  // Imaginary unit (1i since c++17)
  complex<double> I = complex<double>(0., 1.);
  
  // Integration properties
  static const unsigned int NGaussKronrod = 15; // this has to be constant, to be set at compile time
  unsigned int NMaxIter;
  double epsTolerance;
  int NFuncCalls;
};

#endif // _CoulCorrCalc_h_
