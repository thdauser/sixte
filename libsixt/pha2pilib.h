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


   Copyright 2007-2014 Christian Schmid, FAU
*/

#ifndef PHA2PILIB_H
#define PHA2PILIB_H 1

#include "sixt.h"
#include "event.h"
#include "eventfile.h"
#include "rmf.h"
#include "gsl/gsl_rng.h"

#include <unistd.h>

////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////

// PHA2PI /** FITS file containing the correction from PHA to PI. */
typedef struct {

	/** FILE NAME */
	char* pha2pi_filename;

	/** RANDOM NUMBER GENERATOR */
	gsl_rng *randgen;

	/** FILE CONTENT */
	char* rmffile;

	long nrows, ngrades;

	long* pha;
	double** pien;
	double** pilow;
	double** pihigh;

} Pha2Pi;


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty Pha2Pi data
    structure. */
Pha2Pi* getPha2Pi(int* const status);

/** Destructor. */
void freePha2Pi(Pha2Pi** const p2p);

/** Initialize RNG and Load Pha2Pi structure from File. */
Pha2Pi* initPha2Pi(const char* const filename,
		const unsigned int seed,
		int* const status);

/**  Binary search for to find interpolation interval
 *   - return value is the bin [ind,ind+1]
 *   - assume list is sorted ascending */
int binary_search_float(float val, float* arr, long n);


/** Do the pha2pi correction on a single event. */
void pha2pi_correct_event(Event* const evt,
		const Pha2Pi* const p2p,
		const struct RMF* const rmf,
		int* const status);

/** Do the pha2pi correction on a eventfile. */
void pha2pi_correct_eventfile(EventFile* const evtfile,
		const Pha2Pi* const p2p,
		const char* RSPPath,
		const char* RESPfile,
		int* const status);


#endif /* PHA2PILIB_H */

