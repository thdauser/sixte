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

#ifndef PHA2PI_H
#define PHA2PI_H 1


#include "sixt.h"
#include "event.h"
#include "eventfile.h"
#include "parinput.h"
#include "pha2pilib.h"

#define TOOLSUB pha2pi_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters {
  char EvtFile[MAXFILENAME];
  char Pha2Pi[MAXFILENAME];
  char RSPPath[MAXFILENAME];
  char RESPfile[MAXFILENAME];

  int Seed;
  int clobber;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL.
void pha2pi_getpar(struct Parameters* const par, int* const status);



#endif /* pha2pi_H */
