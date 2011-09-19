#ifndef RUNSIXT_H
#define RUNSIXT_H 1

#include "sixt.h"

#include "attitudecatalog.h"
#include "eventlistfile.h"
#include "gendet.h"
#include "impactlistfile.h"
#include "phdet.h"
#include "phgen.h"
#include "phimg.h"
#include "photonlistfile.h"
#include "phproj.h"
#include "sourcecatalog.h"
#include "vector.h"

#define TOOLSUB runsixt_main
#include "headas_main.c"


struct Parameters {
  char OutputStem[MAXFILENAME];
  char PhotonList[MAXFILENAME];
  char ImpactList[MAXFILENAME];
  char EventList[MAXFILENAME];
  char Mission[MAXMSG];
  char Instrument[MAXMSG];
  char Mode[MAXMSG];
  char XMLFile[MAXFILENAME];
  char Attitude[MAXFILENAME];

  /** [deg] */
  float RA, Dec;

  char Simput[MAXFILENAME];

  double MJDREF;
  double TIMEZERO;
  double Exposure;

  int Seed;
  
  char clobber;
};


int runsixt_getpar(struct Parameters* const par);


#endif /* RUNSIXT_H */

