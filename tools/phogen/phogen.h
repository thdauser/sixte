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

#ifndef PHOGEN_H
#define PHOGEN_H 1

#include "sixt.h"

#include "attitude.h"
#include "geninst.h"
#include "phgen.h"
#include "photonfile.h"
#include "simput.h"
#include "vector.h"
#include "simput.h"
#include "sourcecatalog.h"

#define TOOLSUB phogen_main
#include "headas_main.c"


struct Parameters {
  char PhotonList[MAXFILENAME];
  char Mission[MAXMSG];
  char Instrument[MAXMSG];
  char Mode[MAXMSG];
  char XMLFile[MAXFILENAME];
  char Attitude[MAXFILENAME];

  /** [deg] */
  float RA, Dec, rollangle;

  char Simput[MAXFILENAME];

  double MJDREF;
  double TSTART;
  double Exposure;
  double dt;

  int Seed;

  char clobber;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


int phogen_getpar(struct Parameters* parameters);


#endif /* PHOGEN_H */
