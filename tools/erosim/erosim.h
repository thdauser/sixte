#ifndef EROSIM_H
#define EROSIM_H 1

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
#include "phpat.h"
#include "phproj.h"
#include "sourcecatalog.h"
#include "vector.h"

#define TOOLSUB erosim_main
#include "headas_main.c"


/** Maximum number of SIMPUT catalogs. */
#define MAX_N_SIMPUT 6


struct Parameters {
  char Prefix[MAXFILENAME];
  char PhotonList[MAXFILENAME];
  char ImpactList[MAXFILENAME];
  char EventList[MAXFILENAME];
  char PatternList[MAXFILENAME];
  char Attitude[MAXFILENAME];
  char GTIFile[MAXFILENAME];
  char ProgressFile[MAXFILENAME];

  char XMLFile[MAXFILENAME];
  char XMLFile1[MAXFILENAME];
  char XMLFile2[MAXFILENAME];
  char XMLFile3[MAXFILENAME];
  char XMLFile4[MAXFILENAME];
  char XMLFile5[MAXFILENAME];
  char XMLFile6[MAXFILENAME];
  char XMLFile7[MAXFILENAME];

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


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


int erosim_getpar(struct Parameters* const par);


#endif /* EROSIM_H */

