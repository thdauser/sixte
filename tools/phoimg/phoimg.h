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

#ifndef PHOIMG_H
#define PHOIMG_H 1


#include "sixt.h"
#include "attitude.h"
#include "geninst.h"
#include "impactfile.h"
#include "photonfile.h"
#include "phimg.h"

#define TOOLSUB phoimg_main
#include "headas_main.c"


struct Parameters {
  char PhotonList[MAXFILENAME];
  char ImpactList[MAXFILENAME];
  char XMLFile[MAXFILENAME];
  char Attitude[MAXFILENAME];

  /** [deg] */
  float RA, Dec, rollangle;

  double MJDREF;
  double TSTART;
  double Exposure;

  int Seed;

  char clobber;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int phoimg_getpar(struct Parameters* parameters);


#endif /* PHOIMG_H */
