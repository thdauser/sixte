#ifndef COMADET_H
#define COMADET_H 1


#include "sixt.h"
#include "impact.h"
#include "impactfile.h"
#include "comadetector.h"


#define TOOLSUB comadet_main
#include "headas_main.c"


struct Parameters {
  char ImpactList[MAXMSG];
  char EventList[MAXMSG];
  char EventListTemplate[MAXMSG];

  /** Detector width in [pixel]. */
  //has to include all gaps
  int width;
  /** Width of one detector pixel in [m]. */
  double pixelwidth;

  /**length of DCU, gap between 2 DCU's and gap between two DCA's [m]. */
  //only works for 2x2 DCU's separated by DCU_gap, followed by DCA_gap
  double DCU_length;
  double DCU_gap;
  double DCA_gap;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int comadet_getpar(struct Parameters* parameters);


#endif /* COMADET_H */
