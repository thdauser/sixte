#ifndef MAKELC_H
#define MAKELC_H 1

#include "sixt.h"
#include "lad.h"
#include "ladeventlistfile.h"
#include "ladevent.h"

#define TOOLSUB makelc_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char EventList[MAXFILENAME];
  char LightCurve[MAXFILENAME];

  /** Start time [s]. */
  double TSTART;
  
  /** Length of the light curve [s]. */
  double length; 

  /** Time resolution [s]. */
  double dt; 

  /** Lower and upper boundary of the regarded energy band [keV]. */
  float Emin, Emax;

  char clobber;
};


int makelc_getpar(struct Parameters *par);


#endif /* MAKELC_H */

