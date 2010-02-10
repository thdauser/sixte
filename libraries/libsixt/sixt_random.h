#ifndef SIXT_RANDOM_H
#define SIXT_RANDOM_H 1

#include "sixt.h"


/** This routine returns a random number created by the HEAdas random
    number generator. It is basically a wrapper around the HEAdas
    routine HDmtDrand(). The return values lie in the interval
    [0,1). */
inline double sixt_get_random_number();


/** Returns a random value on the basis of an exponential distribution
    with a given average distance. In the simulation this function is
    used to calculate the temporal differences between individual
    photons from a source. The photons have Poisson statistics. */
inline double rndexp(double avg);


/** Determine 2 (!) Gaussian distributed random numbers using the
    Box-Muller method (Gould & Tobochnik, p. 432). The standard
    deviation of the random numbers is 1. */
inline void get_gauss_random_numbers(double* x, double* y);


#endif /* SIXT_RANDOM_H */

