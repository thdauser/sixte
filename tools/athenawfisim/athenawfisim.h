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

#ifndef ATHENAPWFISIM_H
#define ATHENAPWFISIM_H 1

#include "sixt.h"

#include "attitude.h"
#include "eventfile.h"
#include "geninst.h"
#include "gentel.h"
#include "gti.h"
#include "impactfile.h"
#include "phdet.h"
#include "phgen.h"
#include "phimg.h"
#include "photonfile.h"
#include "pha2pilib.h"
#include "phpat.h"
#include "phproj.h"
#include "sourcecatalog.h"
#include "vector.h"

#define TOOLSUB athenapwfisim_main
#include "headas_main.c"


/** Maximum number of SIMPUT catalogs. */
#define MAX_N_SIMPUT 6


struct Parameters {
  char Prefix[MAXFILENAME];
  char PhotonList[MAXFILENAME];
  char ImpactList[MAXFILENAME];
  char EvtFile[MAXFILENAME];
  char RawData[MAXFILENAME];
  char* Attitude;
  char GTIfile[MAXFILENAME];
  char ProgressFile[MAXFILENAME];

  char XMLFile0[MAXFILENAME];
  char XMLFile1[MAXFILENAME];
  char XMLFile2[MAXFILENAME];
  char XMLFile3[MAXFILENAME];

  char Background;

  /** Telescope pointing direction [deg]. */
  float RA, Dec;

  char Simput[MAXFILENAME];
  char Simput2[MAXFILENAME];
  char Simput3[MAXFILENAME];
  char Simput4[MAXFILENAME];
  char Simput5[MAXFILENAME];
  char Simput6[MAXFILENAME];

  double MJDREF;
  double TSTART;
  double Exposure;
  double dt;

  int Seed;

  /** Skip invalid patterns when producing the output file. */
  char SkipInvalids;

  char clobber;
};


int athenapwfisim_getpar(struct Parameters* const par);


#endif /* ATHENAPWFISIM_H */
