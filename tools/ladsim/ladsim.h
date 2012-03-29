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
#include "ladsignallistfile.h"
#include "linkedladsiglist.h"
#include "rmf.h"
#include "sourcecatalog.h"
#include "vector.h"

#define TOOLSUB ladsim_main
#include "headas_main.c"

// Flag whether the Open Area Ratio of the collimator on the LAD
// should be determined at the beginning of the simulation.
#define LAD_OAR 1

struct Parameters {
  char Prefix[MAXFILENAME];
  char PhotonList[MAXFILENAME];
  char ImpactList[MAXFILENAME];
  char SignalList[MAXFILENAME];
  char EventList[MAXFILENAME];
  char XMLFile[MAXFILENAME];
  char Attitude[MAXFILENAME];

  /** [deg] */
  float RA, Dec;

  char Simput[MAXFILENAME];

  double MJDREF;
  double TIMEZERO;
  double Exposure;
  double dt;

  int Seed;
  
  char clobber;
};


int ladsim_getpar(struct Parameters* const par);


#endif /* LADSIM_H */

