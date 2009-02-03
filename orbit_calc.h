//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// general information:
// file:             orbit.h (header file of 'orbit_calc.c')
// project:          eROSITA - NRTA - simulation
// author:           Christian Schmid
// date:             01/2008
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// description:
// The code contains several constants, structures and functions to calculate the orbit 
// of a satellite including perturbations by the particular shape of the earth (see orbit.c). 
// The routines are implemented for maximum speed, so there are some redundant variables.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef ORBIT_CALC_H
#define ORBIT_CALC_H 1

#include <math.h>
#include "vector.h"
#include "earth_constants.h"


// structure with the necessary data for the orbit calculation
// (Keplerian orbital elements)
struct orbit_data {
  long double a;     // semimajor axis
  long double p;     // parameter of the ellipse
  long double e;     // eccentricity
  long double i;     // inclination
  long double Omega; // right ascension of ascending node
  long double omega; // argument of perigee
  long double M;     // mean anomaly
  long double n;     // mean motion

  long double dt;    // width of a timestep
};


// structure which contains the position and velocity of the satellite in its orbit
struct orbit_position {
  struct vector r,v;    // vectors \vec{r}(t) and \vec{v}(t)
  long double t;        // time
};



///////////////////////
// functions/methods //
///////////////////////


// orbit calculation:

// This function calculates the initial orbit data and saves it to the corresponding structure:
void orbit_init(double a, double e, double i, double Omega, double omega, double M, struct orbit_data *odata);


// This function performs a timestep, i.e. it calculates the next position of the satellite
// on its orbit, considering first order perturbation terms, i.e. only J2.
void orbit_step_J2(struct orbit_data *odata, struct orbit_position *oposition);


// This function performs a timestep, i.e. it calculates the next position of the satellite
// on its orbit, considering higher order perturbation terms J2, J2^2, J3, J4.
void orbit_step_J234(struct orbit_data *odata, struct orbit_position *oposition);


// This function performs a timestep, i.e. it calculates the next position of the satellite
// on its orbit, considering higher order perturbation terms J2, J2^2, J3, J4.
// The calculation used follows the approach of Lyddane.
void orbit_step_J234t(struct orbit_data *odata, struct orbit_position *oposition);


// Function solves the Kepler-equation for given mean anomaly 'M' and
// given eccentricity 'e' with the Newton algorithm.
// returns the eccentric anomaly 'E'
long double kepler_equation(const long double M, const long double e);



#endif
