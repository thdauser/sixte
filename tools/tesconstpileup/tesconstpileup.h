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


   Copyright 2014 Philippe Peille, IRAP
*/

#ifndef TESCONSTPILEUP_H
#define TESCONSTPILEUP_H 1

#include <sys/time.h>
#include "sixt.h"
#include "simput.h"
#include "gti.h"
#include "advdet.h"
#include "pixelimpactfile.h"

#define TOOLSUB tesconstpileup_main
#include "headas_main.c"

////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////

struct Parameters {
  char PixImpList[MAXFILENAME];
  char XMLFile[MAXFILENAME];
  
  char telescop[MAXMSG];
  char instrume[MAXMSG];
  char filter[MAXMSG];
  char ancrfile[MAXMSG];
  char respfile[MAXMSG];
  
  double mjdref;
  double timezero;
  double tstop;

  int pulseDistance; // separation in samples between the first 2 pileup events
  int pulseDistance2; // separation in samples between the middle and third events
  double energy; // energy of the primary photons
  double energy2;// energy of the secondary photons
  double energy3;// energy of the tertiary photons
  double offset;// offset in time bin

  //Necessary to ensure exactly 2/3 pulses per trigger
  int preBufferSize;
  int triggerSize;
  
  char clobber;
  char history;
  
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int getpar(struct Parameters* const parameters);

#endif /* TESCONSTPILEUP_H */
