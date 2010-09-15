#ifndef SIXT_RANDOM_H
#define SIXT_RANDOM_H 1

#include "sixt.h"


/** This routine returns a random number created by the HEAdas random
    number generator. It is basically a wrapper around the HEAdas
    routine HDmtDrand(). The return values lie in the interval
    [0,1). The function requires HDmtInit() to be called once before
    usage in order to initialize the HEAdas random number
    generator. When the random number generator is not needed any
    more, it can be realeased with HDmtFree(). Information can be
    found in the HEAdas developer's guide or directly in the source
    files 'headas_rand.h' and 'headas_rand.c'. */
double sixt_get_random_number();


/** Returns a random value on the basis of an exponential distribution
    with a given average distance. In the simulation this function is
    used to calculate the temporal differences between individual
    photons from a source. The photons have Poisson statistics. */
double rndexp(const double avg);


/** Determine 2 (!) Gaussian distributed random numbers using the
    Box-Muller method (Gould & Tobochnik, p. 432). The standard
    deviation of the random numbers is 1. */
void get_gauss_random_numbers(double* const x, double* const y);


#endif /* SIXT_RANDOM_H */

