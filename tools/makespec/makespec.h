#ifndef MAKESPEC_H
#define MAKESPEC_H 1

#include "sixt.h"
#include "rmf.h"
#include "arf.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif

#define TOOLSUB makespec_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char EventList[MAXFILENAME];
  char Spectrum[MAXFILENAME];
  char RSPPath[MAXFILENAME];

  int Seed;
  
  char clobber;
};


int makespec_getpar(struct Parameters *par);


#endif /* MAKESPEC_H */

