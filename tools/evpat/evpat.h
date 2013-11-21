#ifndef EVPAT_H
#define EVPAT_H 1


#include "sixt.h"
#include "event.h"
#include "eventfile.h"
#include "geninst.h"
#include "gti.h"
#include "phpat.h"

#define TOOLSUB evpat_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters {
  char Mission[MAXMSG];
  char Instrument[MAXMSG];
  char Mode[MAXMSG];
  char XMLFile[MAXFILENAME];

  char EventList[MAXFILENAME];
  char PatternList[MAXFILENAME];

  /** Skip invalid patterns when producing the output file. */
  char SkipInvalids;

  char clobber;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL.
int getpar(struct Parameters* const par);


#endif /* EVPAT_H */

