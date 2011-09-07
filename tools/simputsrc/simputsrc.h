#ifndef SIMPUTSRC_H
#define SIMPUTSRC_H 1

#include "sixt.h"
#include "simput.h"

#define TOOLSUB simputsrc_main
#include "headas_main.c"


struct Parameters {

  char Spectrum[MAXFILENAME];
  char Simput[MAXFILENAME];

  /** Reference energy band [keV]. */
  float Emin;
  float Emax;
  float Flux;
  
  char clobber;
};


int simputsrc_getpar(struct Parameters* const par);


#endif /* SIMPUTSRC_H */

