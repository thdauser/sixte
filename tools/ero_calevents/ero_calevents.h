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

#ifndef ERO_CALEVENTS_H
#define ERO_CALEVENTS_H 1


#include "sixt.h"
#include "attitude.h"
#include "event.h"
#include "eventfile.h"
#include "gti.h"
#include "wcs.h"
#include "radec2xylib.h"

#include <libgen.h>

#define TOOLSUB ero_calevents_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters {
  char EvtFile[MAXFILENAME];
  char eroEvtFile[MAXFILENAME];
  char* Attitude;

  int CCDNr;

  /** Projection type (usually SIN). */
  char Projection[MAXMSG];
  /** Right ascension of reference point [deg]. */
  float RefRA;
  /** Declination of reference point [deg]. */
  float RefDec;
  /** Right ascension of pointing [deg]. */
  float RA;
  /** Declination of pointing [deg]. */
  float Dec;
  /** Roll angle of pointing [deg]. */
  float rollangle;

  int usepha;

  char clobber;
};


/** Event entry in an eROSITA calibrated event file. */
typedef struct {
  double time; /* [s] */
  long frame;

  long pha; /* [adu] */
  float pi; /* [eV] */
  float energy; /* [eV] */

  int rawx, rawy;
  int subx, suby;

  double ra, dec; /* [deg] */
  long x, y;

  long flag;

  unsigned int pat_typ;
  unsigned char pat_inf;

  float ev_weight;

  unsigned char ccdnr;

  double recordtime; /* [s] since 1.1.2000 */
  double frametime; /* [s] since 1.1.2000 */
  
} eroCalEvent;


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int getpar(struct Parameters* const parameters);


#endif /* ERO_CALEVENTS_H */
