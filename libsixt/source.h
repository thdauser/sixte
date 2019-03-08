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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef SOURCE_H
#define SOURCE_H 1

#include "sixt.h"

#include "linkedpholist.h"
#include "photon.h"
#include "simput.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Photon energy spectrum of an X-ray source. */
typedef struct {

  /** Coordinates of source position. */
  double ra, dec;

  /** Source extension [rad]. */
  float extension;

  /** Row number of the source in the SimputCtlg. Numbering starts
      at line 1. */
  long row;

  /** Time of the emission of the last photon. */
  double* t_next_photon;

} Source;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. */
Source* newSource(int* const status);

/** Destructor. */
void freeSource(Source** const src);

/** Create photons for a particular source in the specified time
    interval. */
LinkedPhoListElement* getXRayPhotons(Source* const src,
				     SimputCtlg* const simput,
				     const double t0,
				     const double t1,
				     const double mjdref,
				     int* const status);

/** Sort the list of Source objects with the specified number of
    entries with respect to the requested coordinate axis using a
    quick sort algorithm. */
void quicksortSources(Source* const list,
		      const long left,
		      const long right,
		      const int axis);


#endif /* SOURCE_H */
