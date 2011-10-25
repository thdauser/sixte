#ifndef LADLC_H
#define LADLC_H 1

#include "sixt.h"
#include "lad.h"
#include "ladeventlistfile.h"
#include "ladevent.h"

#define TOOLSUB ladlc_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char EventList[MAXFILENAME];
  char LightCurve[MAXFILENAME];
  
  /** Length of the light curve [s]. */
  double length; 

  /** Time resolution [s]. */
  double dt; 
  
  char clobber;
};


int ladlc_getpar(struct Parameters *par);


#endif /* LADLC_H */

