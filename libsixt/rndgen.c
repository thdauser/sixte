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

#include "simput.h"
#include "rndgen.h"

// Use a Pseudo RNG for testing
#include<mt19937ar.h>



/** BOOLEAN: 1 if sixt_init_rng was performed, else 0
 *  sixt_destroy_rng sets this boolean back to 0 */
unsigned int SIXT_RNG_INITIALIZED = 0;

int USE_PSEUDO_RNG = 0;

unsigned int sixt_rng_is_initialized() {
	return SIXT_RNG_INITIALIZED;
}

unsigned int sixt_use_pseudo_rng() {
	return USE_PSEUDO_RNG;
}

double sixt_get_random_number(int* const status){

	// Return a random value out of the interval [0,1).
	if (USE_PSEUDO_RNG==1){
		// use the Pseudo RNG (MT19937) in the [0,1) interval
		return genrand_real2();
	}


#ifdef USE_RCL
  // Use the Remeis random number server running on leo.
  // The library for accessing the server is maintained
  // by Fritz-Walter Schwarm.
  int rcl_status=0;
  double rand=rcl_rand_net_get_double(NULL, NULL, &rcl_status);

  if(RCL_RANDOM_SUCCESS!=rcl_status) {
    SIXT_ERROR("failed getting random number from RCL");
    *status=EXIT_FAILURE;
  }

  return(rand);

#else

  // Use the HEAdas random number generator.
  return(HDmtDrand());

  // Status variable is not needed.
  (void)(*status);
#endif
}


void sixt_init_rng(const unsigned int seed, int* const status) {
  // Initialize HEAdasS random number generator.
  // Note that this has to be done in any case, even
  // if the RCL random number server is used, because
  // the HEAdas routines (like heasp) rely on HDmtDrand().

	if( SIXT_RNG_INITIALIZED == 1 ){
		return;
	}

	setSimputRndGen( &sixt_get_random_number);

	if(getenv("SIXTE_USE_PSEUDO_RNG")!=NULL) {
		USE_PSEUDO_RNG = 1;
		SIXT_WARNING(" using PSEUDO RANDOM NUMBERS (should be only used for testing)!");
		init_genrand(seed);
		SIXT_RNG_INITIALIZED=1;
	}

	HDmtInit(seed);
	SIXT_RNG_INITIALIZED=1;

#ifdef USE_RCL

  // Call the RCL random number generator specifying the
  // server and method.
  int rcl_status=0;
  rcl_rand_net_get_double("draco", "rand", &rcl_status);

  if(RCL_RANDOM_SUCCESS!=rcl_status) {
    SIXT_ERROR("failed getting random number from RCL");
    *status=EXIT_FAILURE;
  }

#else

  // The status variable is not used.
  (void)(*status);

#endif
}


void sixt_destroy_rng()
{
  // Release HEADAS random number generator:
	if( SIXT_RNG_INITIALIZED==1 ){
		HDmtFree();
		SIXT_RNG_INITIALIZED=0;
	}

  // remove preference to use PSEUDO RNG
	USE_PSEUDO_RNG = 0;
}


void sixt_get_gauss_random_numbers(double* const x,
				   double* const y,
				   int* const status)
{
  double sqrt_2rho=sqrt(-log(sixt_get_random_number(status))*2.);
  CHECK_STATUS_VOID(*status);
  double phi=sixt_get_random_number(status)*2.*M_PI;
  CHECK_STATUS_VOID(*status);

  *x=sqrt_2rho * cos(phi);
  *y=sqrt_2rho * sin(phi);
}


double rndexp(const double avgdist, int* const status)
{
  assert(avgdist>0.);

  double rand=sixt_get_random_number(status);
  CHECK_STATUS_RET(*status, 0.);
  if (rand<1.e-15) {
    rand=1.e-15;
  }

  return(-log(rand)*avgdist);
}

