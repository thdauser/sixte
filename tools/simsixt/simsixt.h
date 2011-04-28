#ifndef SIMSIXT_H
#define SIMSIXT_H 1

#include "sixt.h"

#include "attitudecatalog.h"
#include "gendet.h"
#include "impactlistfile.h"
#include "phdet.h"
#include "phgen.h"
#include "phimg.h"
#include "photonlistfile.h"
#include "vector.h"
#include "sourcecatalog.h"

#define TOOLSUB simsixt_main
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

  char data_path[MAXFILENAME];
};


int simsixt_getpar(struct Parameters* const par);


#endif /* SIMSIXT_H */

