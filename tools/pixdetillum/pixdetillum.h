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


   Copyright 2014 Thorsten Brand, FAU
*/

#ifndef DETPIXILLUM_H
#define DETPIXILLUM_H 1

#include <sys/time.h>
#include "sixt.h"
#include "simput.h"
#include "gti.h"
#include "advdet.h"
#include "pixelimpactfile.h"

#define TOOLSUB pixdetillum_main
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
  double tstart;
  double tstop;

  int Nactive;
  int nlo;
  int nhi;
  
  double rate;
  double energy;
  
  char clobber;
  char history;
  
  unsigned long int seed;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int getpar(struct Parameters* const parameters);

#endif /* DETPIXILLUM_H */