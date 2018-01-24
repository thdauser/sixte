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

#ifndef XIFUPIPELINE_H
#define XIFUPIPELINE_H 1

#include "sixt.h"

#include "attitude.h"
#include "eventfile.h"
#include "geninst.h"
#include "gti.h"
#include "impactfile.h"
#include "phdet.h"
#include "phgen.h"
#include "phimg.h"
#include "photonfile.h"
#include "phpat.h"
#include "phproj.h"
#include "sourcecatalog.h"
#include "vector.h"
#include "advdet.h"
#include "pixelimpactfile.h"
#include "tesinitialization.h"
#include "testrigger.h"
#include "teseventlist.h"
#include "crosstalk.h"
#include "parinput.h"
#include "grading.h"

#define TOOLSUB xifupipeline_main
#include "headas_main.c"


/** Maximum number of SIMPUT catalogs. */
#define MAX_N_SIMPUT 6


struct Parameters {
  char Prefix[MAXFILENAME];
  char PhotonList[MAXFILENAME];
  char ImpactList[MAXFILENAME];
  char PixImpactList[MAXFILENAME];
  char TesTriggerFile[MAXFILENAME];
  char EvtFile[MAXFILENAME];
  char XMLFile[MAXFILENAME];
  char AdvXml[MAXFILENAME];
  char* Attitude;
  char GTIfile[MAXFILENAME];
  char ProgressFile[MAXFILENAME];
  char PulseTemplateFile[MAXFILENAME];
  char OptimalFilterFile[MAXFILENAME];

  char Mission[MAXMSG];
  char Instrument[MAXMSG];
  char Mode[MAXMSG];

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

  /** Skip invalid patterns when producing the output pattern file. */
  char SkipInvalids;

  /** Tes simulation parameters*/
  int triggerSize;
  int preBufferSize;
  int PulseLength;
  int EventListSize;
  int NormalExclusion;
  int DerivateExclusion;

  double tstop;
  double Threshold;
  double Calfac;
  double SaturationValue;

  char writeStreamFile;
  char Reconstruct;
  char WriteRecordFile;
  char Identify;
  char UseRMF;
  char ProjCenter;

  int doCrosstalk;
  int saveCrosstalk;
  float elec_ctk_scaling;

  /** TDM related constants*/
  int tdm;

  char history;
  char clobber;


};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


int xifupipeline_getpar(struct Parameters* const par);

/** Copies the parameters contained in the local parameter structure into the
    general TES parameters structure */
void copyParams2GeneralStruct(const struct Parameters partmp, TESGeneralParameters* const par,double tstart,double tstop);

/** Iterates over the piximpact file to see which pixels were actually hit */
void getListPixelsHit(PixImpFile* pixilf,int** list_pixels,int npix,int* const status);

#endif /* XIFUPIPELINE_H */

