#ifndef SIMPUTSRC_H
#define SIMPUTSRC_H 1

#include "sixt.h"
#include "simput.h"

#define TOOLSUB simputsrc_main
#include "headas_main.c"


struct Parameters {

  /** Power law. */
  float plPhoIndex;
  float plFlux;

  /** Black body temperature [keV]. */
  float bbkT;
  float bbFlux;

  /** Line dispersion [keV]. */
  float flSigma;
  float flFlux;

  float rflSpin;
  float rflFlux;

  /** Absorption column [10^22 atoms/cm^2] */
  float nH;

  /** Reference energy band [keV]. */
  float Emin;
  float Emax;

  /** File name of output SIMPUT file. */
  char Simput[MAXFILENAME];
  
  char clobber;
};


int simputsrc_getpar(struct Parameters* const par);


#endif /* SIMPUTSRC_H */

