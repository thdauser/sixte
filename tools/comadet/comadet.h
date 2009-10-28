#ifndef COMADET_H
#define COMADET_H 1


#include "sixt.h"
#include "impact.h"
#include "impactlistfile.h"
#include "comadetector.h"


#define TOOLSUB comadet_main
#include "headas_main.c"


struct Parameters {
  char impactlist_filename[MAXMSG];
  char eventlist_filename[MAXMSG];
  char eventlist_template[MAXMSG];

  /** Detector width in [pixel]. */
  int width;
  /** Width of one detector pixel in [m]. */
  double pixelwidth;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int comadet_getpar(struct Parameters* parameters);


#endif /* COMADET_H */
