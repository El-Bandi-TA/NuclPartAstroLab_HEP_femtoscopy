#include "functions.h"
#define NMIN_SIN 50 // at least this many points in a quarter period of the sine function
#define NMIN_EXP 40 // at least this many points in one 1/e change of the exponential
#define LIMIT 30    // if the exponent in the exponential in the integrand is less than minus LIMIT, stop iteration
#define EPSILON 1e-12 // this is for the terms in the series of the hypergeometric function.
#define EPSILON_ABC 1e-12 // this is for the terms in the series of the hypergeometric function.
#define EPSILON_ASYM 1e-6 // this is for the terms in the asymptotic series ofdthe hypergeometric function.


// for the Lanczos approximation of the gamma function
const double lanczos_coeff[9] = {
0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
const double lanczos_g = 7.;

double factorial(const int n)
{
  if(n < 0)
    return 1e99;
  double factor = 1.;
  double result = 1.;
  for(int i = 0; i < n; i++)
  {
    result *= factor; 
    factor += 1.;
  }
  return result;
}

complex<double> Gamma(const complex<double> z)  //Utilizing the "Lanczos approximation"
{
  if(real(z) < 0.5) return M_PI / sin(M_PI * z) / Gamma( 1. - z);
  complex<double> A_g = lanczos_coeff[0] + lanczos_coeff[1] / z + lanczos_coeff[2] / (z + 1.) + lanczos_coeff[3] / (z + 2.) + lanczos_coeff[4] / (z + 3.) + 
                        lanczos_coeff[5] / (z + 4.) + lanczos_coeff[6] / (z + 5.) + lanczos_coeff[7] / (z + 6.) + lanczos_coeff[8] / (z + 7.);
  return sqrt(2. * M_PI) * pow(z + lanczos_g - 0.5, z - 0.5 ) * exp(0.5 - z - lanczos_g) * A_g;
}

complex<double> digamma(const complex<double> z)  //Utilizing the "Lanczos approximation"
{
  if(real(z) < 0.5) return - M_PI * cos(M_PI * z) / sin(M_PI * z) + digamma( 1. - z);
  complex<double> Avg = lanczos_coeff[1] / SQR(z) + lanczos_coeff[2] / SQR(z+1.) + lanczos_coeff[3] / SQR(z+2.) + lanczos_coeff[4] / SQR(z+3.) + 
                        lanczos_coeff[5] / SQR(z+4.) + lanczos_coeff[6] / SQR(z+5.) + lanczos_coeff[7] / SQR(z+6.) + lanczos_coeff[8] / SQR(z+7.);
  complex<double> A_g = lanczos_coeff[0] + lanczos_coeff[1] / z + lanczos_coeff[2] / (z+1.) + lanczos_coeff[3] / (z+2.) + lanczos_coeff[4] / (z+3.) + 
                        lanczos_coeff[5] / (z+4.) + lanczos_coeff[6] / (z+5.) + lanczos_coeff[7] / (z+6.) + lanczos_coeff[8] / (z+7.);
  return log( z + lanczos_g - 0.5) - lanczos_g  / (z + lanczos_g - 0.5) - Avg / A_g;
}

double Gamma(const double x)
{
  if(x < 0.5) return M_PI / sin(M_PI * x) / Gamma( 1. - x);
  double A_g = lanczos_coeff[0] + lanczos_coeff[1] / x + lanczos_coeff[2] / (x + 1.) + lanczos_coeff[3] / (x + 2.) + lanczos_coeff[4] / (x + 3.) + 
                        lanczos_coeff[5] / (x + 4.) + lanczos_coeff[6] / (x + 5.) + lanczos_coeff[7] / (x + 6.) + lanczos_coeff[8] / (x + 7.);
  return sqrt(2. * M_PI) * pow(x + lanczos_g - 0.5, x - 0.5 ) * exp(0.5 - x - lanczos_g) * A_g;
}

complex<double> Gauss_2F1_series0(const complex<double> a, const complex<double> b, const complex<double> c, const complex<double> z)
{
  complex<double> result(-9999., 0);
  if(abs(z) < 0.9999)
  {
    complex<double> n(0., 0.);
    complex<double> term(1., 0);
    result = term;
    while(abs(term) > EPSILON)
    {
      term = term * z / (n + 1.) * (a + n) * (b + n) / (c + n);
      result += term;
      n += 1. ;
    }
  }
  return result;
}

complex<double> Gauss_2F1_1mZ_noint(const complex<double> a, const complex<double> b, const complex<double> c, const complex<double> z)
{
  complex<double> result(-9999., 0);
  if(abs(1. - z) < 0.9999)
  {
    complex<double> sinabc = sin(M_PI * (c-a-b));
    if(abs(sinabc) > EPSILON_ABC)
      result = M_PI / sinabc * Gamma(c) * ( Gauss_2F1_series0(a, b, a+b+1.-c, 1.-z) / Gamma(a+b+1.-c) / Gamma(c-a) / Gamma(c-b)
                                           - Gauss_2F1_series0(c-a, c-b, c+1.-a-b, 1.-z) / Gamma(a) / Gamma(b) / Gamma(c+1.-a-b) / pow(1.-z, a+b-c));
  }
  return result;
}