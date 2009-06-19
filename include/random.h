/** 
 * Header file of random.c
 */

#ifndef RANDOM_H
#define RANDOM_H 1

#include <stdlib.h>
#include <limits.h>
#include <math.h>

#include "headas.h"
#include "headas_rand.h"


// This routines Creates random numbers using the HEAdas random number generator.
inline double get_random_number();

// Returns a random value on the basis of an exponential distribution with a given average.
// Here this function is used to calculate the temporal differences between individual photons 
// from a source.
inline double rndexp(double avg);



#endif /* RANDOM_H */

