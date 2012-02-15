#ifndef WFMPHOVIGN_H
#define WFMPHOVIGN_H 1


#include "sixt.h"
#include "attitudecatalog.h"
#include "photonlistfile.h"

#define TOOLSUB wfmphovign_main
#include "headas_main.c"


struct Parameters {
  char InputList[MAXFILENAME];
  char OutputList[MAXFILENAME];
  char Attitude[MAXFILENAME];

  /** [deg] */
  float RA, Dec;

  int Seed;
  
  char clobber;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int wfmphovign_getpar(struct Parameters* parameters);


#endif /* WFMPHOVIGN_H */
