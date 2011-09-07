#ifndef SIMPUTSPEC_H
#define SIMPUTSPEC_H 1

#include "sixt.h"
#include "simput.h"

#define TOOLSUB simputspec_main
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

  char Outfile[MAXFILENAME];
  
  char clobber;
};


int simputspec_getpar(struct Parameters* const par);


#endif /* SIMPUTSPEC_H */

