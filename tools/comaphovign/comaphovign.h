#ifndef COMAPHOVIGN_H
#define COMAPHOVIGN_H 1


#include "sixt.h"
#include "attitudecatalog.h"
#include "photonlistfile.h"

#define TOOLSUB comaphovign_main
#include "headas_main.c"


struct Parameters {
  char InputList[MAXFILENAME];
  char OutputList[MAXFILENAME];
  char Attitude[MAXFILENAME];

  /** [deg] */
  float RA, Dec;

  /** Flag whether the time information of the photons should be
      copied to the output file. The time column might be omitted in
      order to save memory. */
  char TimeColumn;

  int Seed;
  
  char clobber;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int comaphovign_getpar(struct Parameters* parameters);


#endif /* COMAPHOVIGN_H */
