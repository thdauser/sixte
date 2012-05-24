#ifndef RUNSIXT_H
#define RUNSIXT_H 1

#include "sixt.h"

#include "attitudecatalog.h"
#include "eventlistfile.h"
#include "gendet.h"
#include "impactlistfile.h"
#include "patternfile.h"
#include "phdet.h"
#include "phgen.h"
#include "phimg.h"
#include "photonlistfile.h"
#include "phpat.h"
#include "phproj.h"
#include "sourcecatalog.h"
#include "vector.h"

#define TOOLSUB runsixt_main
#include "headas_main.c"


/** Maximum number of SIMPUT catalogs. */
#define MAX_N_SIMPUT 5


struct Parameters {
  char Prefix[MAXFILENAME];
  char PhotonList[MAXFILENAME];
  char ImpactList[MAXFILENAME];
  char EventList[MAXFILENAME];
  char PatternList[MAXFILENAME];
  char XMLFile[MAXFILENAME];
  char Attitude[MAXFILENAME];
  char ProgressFile[MAXFILENAME];

  char Mission[MAXMSG];
  char Instrument[MAXMSG];
  char Mode[MAXMSG];

  /** [deg] */
  float RA, Dec;

  char Simput[MAXFILENAME];
  char Simput2[MAXFILENAME];
  char Simput3[MAXFILENAME];
  char Simput4[MAXFILENAME];
  char Simput5[MAXFILENAME];

  double MJDREF;
  double TIMEZERO;
  double Exposure;
  double dt;

  int Seed;
  
  char clobber;
};


int runsixt_getpar(struct Parameters* const par);


#endif /* RUNSIXT_H */

