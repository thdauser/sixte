#ifndef TES_GRADES_H
#define TES_GRADES_H 1

#include "sixt.h"
#include "eventfile.h"


#define TOOLSUB tes_grades_main
#include "headas_main.c"


struct Parameters{
  char EventList[MAXFILENAME];

  /** Characteristic time unit of the TES microcalorimeter. */
  double TimeUnit;
  int PreTrigger;
  int PostTrigger;
  double PileupTime;

  char clobber;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


int tes_grades_getpar(struct Parameters* par);


#endif /* TES_GRADES_H */
