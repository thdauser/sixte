/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2019 Thomas Dauser, FAU
*/

#ifndef RNDGEN_H
#define RNDGEN_H 1

#include "headas_rand.h"
#include "simput.h"
#include "sixt.h"

#ifdef USE_RCL
// Use the Remeis random number server.
#include <rcl.h>
#endif


/** Return value of SIXT_RNG_INITIALIZED. */
unsigned int sixt_rng_is_initialized();

/** Return value of USE_PSEUDO_RNG. */
unsigned int sixt_use_pseudo_rng();


/** This routine returns a random number. The values are either
    obtained from the Remeis random number server or are created by
    the HEAdas random number generator. The routine is basically a
    wrapper around the respective library routines, either the
    rcl_rand_ndg() or the HEAdas routine HDmtDrand(). The return value
    lies in the interval [0,1).

    In case the HEAdas random number generator is used, the function
    requires HDmtInit() to be called once before usage for
    initialization. When the HEAdas random number generator is not
    needed any more, it can be realeased with HDmtFree(). Information
    can be found in the HEAdas developer's guide or directly in the
    source files 'headas_rand.h' and 'headas_rand.c'. */
double sixt_get_random_number(int* const status);

/** Initialize the random number generator. */
void sixt_init_rng(const unsigned int seed, int* const status);

/** Clean up the random number generator. */
void sixt_destroy_rng();

/** This routine produces two Gaussian distributed random numbers. The
    standard deviation of the Gaussian distribution sigma is assumed
    to be unity. The two numbers are returned via the pointer function
    arguments. */
void sixt_get_gauss_random_numbers(double* const x,
				   double* const y,
				   int* const status);

/** Returns a random value on the basis of an exponential distribution
    with a given average distance. In the simulation this function is
    used to calculate the temporal differences between individual
    photons from a source. The photons have Poisson statistics. */
double rndexp(const double avg, int* const status);


#endif /* RNDGEN_H */
