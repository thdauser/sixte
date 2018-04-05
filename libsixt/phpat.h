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

#ifndef PHPAT_H
#define PHPAT_H 1

#include "sixt.h"
#include "event.h"
#include "eventfile.h"
#include "gendet.h"
#include "gsl/gsl_rng.h"


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


void phpat(GenDet* const det,
	   const EventFile* const src,
	   EventFile* const dest,
		 const char* picorr_file,
	   const char skip_invalids,
	   int* const status);

// PHA2PI /** FITS file containing the correction from PHA to PI. */
typedef struct {

  /** Total number of rows in the file. */
  long nrows, ngrades;

  double* pha;
  double** pien;
  double** pilow;
  double** pihigh;

} Pha2Pi;

/** Constructor. Returns a pointer to an empty Pha2Pi data
    structure. */
Pha2Pi* getPha2Pi(int* const status);

void printPha2Pi(const Pha2Pi* const p2p);

/** Destructor. */
void freePha2Pi(Pha2Pi** const p2p);

/** Load Pha2Pi structure from File. */
Pha2Pi* loadPha2Pi(const char* const filename,
			 int* const status);

/** Do the pha2pi correction. */
void pha2picorrect(Event* const evt, const Pha2Pi* const p2p, int* const status);

#endif /* PHPAT_H */
