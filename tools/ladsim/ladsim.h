#ifndef LADSIM_H
#define LADSIM_H 1

#include "sixt.h"

#include "attitudecatalog.h"
#include "genericdetector.h"
#include "phdet.h"
#include "phgen.h"
#include "phimg.h"
#include "photonlistfile.h"
#include "phproj.h"
#include "lad.h"
#include "ladeventlistfile.h"
#include "ladimpactlistfile.h"
#include "sourcecatalog.h"
#include "vector.h"

#define TOOLSUB ladsim_main
#include "headas_main.c"


struct Parameters {
  char OutputStem[MAXFILENAME];
  char PhotonList[MAXFILENAME];
  char ImpactList[MAXFILENAME];
  char EventList[MAXFILENAME];
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


int ladsim_getpar(struct Parameters* const par);


#endif /* LADSIM_H */

