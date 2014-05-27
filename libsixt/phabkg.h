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

#ifndef PHABKG_H 
#define PHABKG_H 1

#include "sixt.h"
#include "vignetting.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic square pixel grid. */
typedef struct {
  /** Number of spectral bins / channels. */
  long nbins;

  /** PHA channel numbers. */
  long* channel;

  /** Background event rate distribution [counts/s/bin/m^2]. */
  float* distribution;

  /** Time of the next background event. */
  double tnext;

  /** Telescope vignetting function. */
  Vignetting** vignetting;

  /** Telescope focal length. */
  float* focal_length;

} PHABkg;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Loads the PHA spectrum required for the model from
    the specified file. The PHA file must contain two columns: CHANNEL
    and RATE. The entries in the CHANNEL column correspond to the
    channels in the EBOUNDS extension of the RMF. The entries in the
    RATE column represent the background event rate in this particular
    energy channel per second and per illuminated detector area
    [counts/s/bin/m^2]. */
PHABkg* newPHABkg(const char* const filename, int* const status);

/** Destructor. */
void destroyPHABkg(PHABkg** const phabkg);

/** Determine the PHA value and time of a background event according
    to the specified spectral distribution. The function returns an
    individual event. The time differences between the events are
    exponentially distributed. The average rate is determined by the
    distribution given in the PHA data set multiplied with the
    scaling factor. */
int getPHABkgEvent(PHABkg* const phabkg,
		   /** Scaling factor for the count rate
		       distribution. Must be given in [m^2]. */
		   const float scaling,
		   const double tstart,
		   /** Upper limit for the time of the background
		       event [s]. */
		   const double tstop,
		   double* const t,
		   long* const pha,
		   int* const status);


#endif /* PHABKG_H */
