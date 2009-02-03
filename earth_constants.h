#ifndef EARTH_CONSTANTS_H
#define EARTH_CONSTANTS_H

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

#endif
