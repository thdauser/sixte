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


   Copyright 2007-2014 Christian Schmid, Mirjam Oertel, FAU
*/

#ifndef COMADET_H
#define COMADET_H 1


#include "sixt.h"
#include "impact.h"
#include "impactfile.h"
#include "comadetector.h"


#define TOOLSUB comadet_main
#include "headas_main.c"


struct Parameters {
  char ImpactList[MAXMSG];
  char EventList[MAXMSG];
  char EventListTemplate[MAXMSG];

  //protoMirax-flag: 1=yes, 0=no
  int protoMirax;
  /** Detector width in [pixel]. */
  //has to include all gaps
  int width;
  /** Width of one detector pixel in [m]. */
  double pixelwidth;

  /**length of DCU, gap between 2 DCU's and gap between two DCA's [m]. */
  //only works for 2x2 DCU's separated by DCU_gap, followed by DCA_gap
  double DCU_length;
  double DCU_gap;
  double DCA_gap;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int comadet_getpar(struct Parameters* parameters);


#endif /* COMADET_H */
