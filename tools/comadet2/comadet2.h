#ifndef COMADET_H
#define COMADET_H 1


#include "sixt.h"
#include "impact.h"
#include "impactlistfile.h"
#include "comadetector.h"


#define TOOLSUB comadet2_main
#include "headas_main.c"


struct Parameters {
  char ImpactList[MAXMSG];
  char EventList[MAXMSG];
  char EventListTemplate[MAXMSG];

  /** Detector width in [pixel]. */
  int width;
  /** Width of one detector pixel in [m]. */
  double pixelwidth;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int comadet_getpar2(struct Parameters* parameters);


#endif /* COMADET_H */
