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

#ifndef MAKELC_H
#define MAKELC_H 1

#include "sixt.h"

#define TOOLSUB makelc_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char EventList[MAXFILENAME];
  char LightCurve[MAXFILENAME];

  /** Start time [s]. */
  double TSTART;
  
  /** Length of the light curve [s]. */
  double length; 

  /** Time resolution [s]. */
  double dt; 

  /** Lower and upper boundary of the regarded energy band [keV]. */
  float Emin, Emax;

  /** Lower and upper boundary of the regarded channel range [adu]. */
  long Chanmin, Chanmax;

  char clobber;
};


int makelc_getpar(struct Parameters *par);


#endif /* MAKELC_H */

