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

#ifndef EXPOSURE_MAP_H
#define EXPOSURE_MAP_H 1

#include "sixt.h"

#include "attitude.h"
#include "geninst.h"
#include "phimg.h"
#include "phproj.h"
#include "telescope.h"
#include "vector.h"
#include "check_fov.h"
#include "vignetting.h"
#include "parinput.h"
#include "sys/stat.h"

#define TOOLSUB exposure_map_main
#include "headas_main.c"


////////////////////////////////////////////////////

/* Program parameters */
struct Parameters {
  char *Attitude;    // filename of the attitude file
  char *Vignetting;  // filename of the vignetting file
  char *Exposuremap; // output: exposure map
  char *RawExposuremap; // output: raw exposure map (no vignetting)
  char *ProgressFile;
  char *XMLFile;

  /** Telescope Pointing direction [deg]. */
  float RA, Dec;

  double TSTART;
  double timespan;
  /** Step width for the exposure map calculation [s]. */
  double dt;

  int seed;

  /** Right ascension range [rad]. */
  double ra1 , ra2;
  /** Declination range [rad]. */
  double dec1, dec2;
  /** Number of pixels in right ascension and declination. */
  int ra_bins, dec_bins;

  float fov_diameter;

  /** Projection method (1: AIT, 2: SIN). */
  int projection;

  /** Number of interim maps to be stored. */
  int intermaps;

  int clobber;
};

typedef struct {
	char **xmlarray;
	int n;
}xmlarray;

////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////

int exposure_map_getpar(struct Parameters *parameters);


#endif /* EXPOSURE_MAP_H */

