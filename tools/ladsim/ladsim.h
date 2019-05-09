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

#ifndef LADSIM_H
#define LADSIM_H 1

#include "sixt.h"

#include "attitude.h"
#include "genericdetector.h"
#include "gti.h"
#include "phdet.h"
#include "phgen.h"
#include "phimg.h"
#include "photonfile.h"
#include "phproj.h"
#include "lad.h"
#include "ladeventfile.h"
#include "ladimpactfile.h"
#include "ladsignalfile.h"
#include "ladsignallist.h"
#include "rmf.h"
#include "sourcecatalog.h"
#include "vector.h"

#define TOOLSUB ladsim_main
#include "headas_main.c"


// Flag whether the Open Area Ratio of the collimator on the LAD
// should be determined at the beginning of the simulation.
//#define LAD_OAR 1


struct Parameters {
  char Prefix[MAXFILENAME];
  char PhotonList[MAXFILENAME];
  char ImpactList[MAXFILENAME];
  char SignalList[MAXFILENAME];
  char EvtFile[MAXFILENAME];
  char XMLFile[MAXFILENAME];
  char Attitude[MAXFILENAME];
  char Simput[MAXFILENAME];
  char ProgressFile[MAXFILENAME];

  /** [deg] */
  float RA, Dec, rollangle;

  double MJDREF;
  double TSTART;
  double Exposure;
  double dt;

  int Seed;

  char clobber;
};


int ladsim_getpar(struct Parameters* const par);


#endif /* LADSIM_H */
