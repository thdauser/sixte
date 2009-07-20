#ifndef SIXT_RANDOM_H
#define SIXT_RANDOM_H 1

#include "sixt.h"


/** This routine creates random numbers using the HEAdas random number generator.
 * The return values lie in the interval [0,1). */
inline double get_random_number();


/** Returns a random value on the basis of an exponential distribution 
 * with a given average distance. In the simulation this function is used to calculate 
 * the temporal differences between individual photons from a source. The photons have 
 * Poisson statistics. */
inline double rndexp(double avg);


#endif /* SIXT_RANDOM_H */

