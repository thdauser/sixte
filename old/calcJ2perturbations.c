#include <stdio.h>
#include <math.h>

// earth constants:
const long double mu = 398601.3;   // G*Msolar; units: (km^3/s^2)
const long double R_e = 6378.140;  // radius of the earth (mean equatorial radius)
  // physical constants for the earth (Flury, p. 69):
  //      mu  = 398601.3 km^3/s^2
  //      M_e = 5.976 * 10^27 g
  //      R_e = 6378.140 km  (mean equatorial radius)

// Jeffrey coefficients (see Flury p. 70):
const long double J2 = 1082.6268E-6;
const long double J3 =   -2.5356E-6;
const long double J4 =   -1.6234E-6;
//const double J5 =   -0.2276E-6;
//const double J6 =    0.5434E-6;

//////////////////////////////////////////////////////////////////////////////////////////////////////
// This program calculates the J2 orbit perturbations analytically for all values of i \in [0,180]. //
//////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  const double a = 6958.140;
  const double e = 0.0;
  const double timespan = 24. * 3600.;

  int i;
  double inc;
  double factor;

  // calculate prefactor
  factor = 1.5 * sqrt(mu/pow(a,3.)) * J2 /pow(1-pow(e,2.),2.) * pow(R_e/a,2.) * timespan * 180./M_PI;

  for (i=0;i<=180;i++) {
    // calculate the inclination in radians
    inc = (double)i * M_PI/180.;
    // print the inclination to STDOUT
    printf("%lf\t", inc*180./M_PI);
    // \Delta \Omega:
    printf("%lf\t", - factor * cos(inc));
    // \Delta \omega:
    printf("%lf\t", 0.5 * factor * (5. * pow(cos(inc),2.)-1.));
    // \Delta M:
    printf("%lf\n", 0.5 * factor * (3. * pow(cos(inc),2.)-1.) + sqrt(mu/pow(a,3.)) * timespan * 180./M_PI);
  }

  return(0);
}
