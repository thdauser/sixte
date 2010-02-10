#include "sixt_random.h"


inline double sixt_get_random_number()
{
  // Annotation:
  // The function HDmtRand() returns a value of the type 
  // "unsigned long int".

  // Return a value out of the interval [0,1):
  //  return(((double)HDmtRand())/ULONG_MAX);

  return(HDmtDrand());
}



inline double rndexp(double avgdist)
{
  double rand = sixt_get_random_number();
  if (rand < 1.E-15) {
    rand = 1.E-15;
  }
  
  return(-log(rand)*avgdist);
}



inline void get_gauss_random_numbers(double* x, double* y)
{
  double sqrt_2rho = sqrt(-log(sixt_get_random_number())*2.);
  double phi = sixt_get_random_number()*2.*M_PI;

  *x = sqrt_2rho * cos(phi);
  *y = sqrt_2rho * sin(phi);
}


