#ifndef ANALYSE_XMS_EVENTS_H
#define ANALYSE_XMS_EVENTS_H 1

#include "sixt.h"
#include "eventlistfile.h"
#include "patternfile.h"


#define TOOLSUB analyse_xms_events_main
#include "headas_main.c"


struct Parameters{
  char EventList[MAXFILENAME];
  char PatternList[MAXFILENAME];

  /** Characteristic time unit of the TES microcalorimeter. */
  double TimeUnit;
  int PreTrigger;
  int PostTrigger;
  double PileupTime;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


int analyse_xms_events_getpar(struct Parameters* par);


#endif /* ANALYSE_XMS_EVENTS_H */
