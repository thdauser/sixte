#include "sixt_random.h"


double sixt_get_random_number()
{
  // Return a value out of the interval [0,1):
  return(HDmtDrand());
}



double rndexp(const double avgdist)
{
  assert(avgdist>0.);

  double rand = sixt_get_random_number();
  if (rand < 1.E-15) {
    rand = 1.E-15;
  }

  return(-log(rand)*avgdist);
}



void get_gauss_random_numbers(double* const x, double* const y)
{
  double sqrt_2rho = sqrt(-log(sixt_get_random_number())*2.);
  double phi = sixt_get_random_number()*2.*M_PI;

  *x = sqrt_2rho * cos(phi);
  *y = sqrt_2rho * sin(phi);
}


