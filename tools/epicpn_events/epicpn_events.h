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


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

#ifndef EPICPN_EVENTS_H
#define EPICPN_EVENTS_H 1


#include "sixt.h"
#include "attitude.h"
#include "event.h"
#include "eventfile.h"
#include "wcs.h"

#define TOOLSUB epicpn_events_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters {
  char EvtFile[MAXFILENAME];
  char EPICpnEventList[MAXFILENAME];

  char clobber;
};


/** Event entry in an EPIC-pn event file. */
typedef struct {
  /* [s] */
  double time;

  /* [pixel] */
  int rawx, rawy;
  /* [0.05arcsec] */
  int detx, dety;
  /* [0.05arcsec] */
  long x, y;

  /* [adu] */
  int pha;
  /* [eV] */
  int pi;

  long flag;
  char pattern;

} EPICpnEvent;


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int getpar(struct Parameters* const parameters);


#endif /* EPICPN_EVENTS_H */
