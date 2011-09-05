#ifndef SIMPUTSPEC_H
#define SIMPUTSPEC_H 1

#include "sixt.h"
#include "simput.h"

#define TOOLSUB simputspec_main
#include "headas_main.c"


struct Parameters {
  float plIndex;
  float plNorm;
  float bbkT;
  float bbNorm;
  float flSigma;
  float flNorm;
  float rflSpin;
  float rflNorm;

  char Outfile[MAXFILENAME];
  
  char clobber;
};


int simputspec_getpar(struct Parameters* const par);


#endif /* SIMPUTSPEC_H */

