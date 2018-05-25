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

#ifndef MAKESPEC_H
#define MAKESPEC_H 1

#include "sixt.h"
#include "rmf.h"
#include "arf.h"
#include "gti.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif

#define TOOLSUB makespec_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char EvtFile[MAXFILENAME];
  char Spectrum[MAXFILENAME];
  char RSPPath[MAXFILENAME];
  char EventFilter[MAXFILENAME];  
  char ARFfile[MAXFILENAME];
  char RMFfile[MAXFILENAME];

  int Seed;
  
  char clobber;
  char usepha;
};


int makespec_getpar(struct Parameters *par);


#endif /* MAKESPEC_H */

