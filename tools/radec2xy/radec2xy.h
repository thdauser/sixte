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

#ifndef RADEC2XY_H
#define RADEC2XY_H 1


#include "sixt.h"
#include "event.h"
#include "eventfile.h"
#include "wcs.h"
#include "parinput.h"
#include "radec2xylib.h"

#define TOOLSUB radec2xy_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters {
  char EvtFile[MAXFILENAME];
  /** Projection type (usually SIN). */
  char Projection[MAXFILENAME];
  /** Right ascension of reference point [deg]. */
  float RefRA;
  /** Declination of reference point [deg]. */
  float RefDec;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
void radec2xy_getpar(struct Parameters* const parameters, int* const status);


#endif /* ERO_CALEVENTS_H */
