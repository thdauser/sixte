#ifndef PHOVIGN_H
#define PHOVIGN_H 1


#include "sixt.h"
#include "attitudecatalog.h"
#include "photonlistfile.h"

#define TOOLSUB phovign_main
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
int phovign_getpar(struct Parameters* parameters);


#endif /* PHOVIGN_H */
