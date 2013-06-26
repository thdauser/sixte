#ifndef ERO_CALEVENTS_H
#define ERO_CALEVENTS_H 1


#include "sixt.h"
#include "attitude.h"
#include "gti.h"
#include "pattern.h"
#include "patternfile.h"
#include "wcs.h"

#define TOOLSUB ero_calevents_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters {
  char PatternList[MAXFILENAME];
  char eroEventList[MAXFILENAME];
  char GTIFile[MAXFILENAME];
  char Attitude[MAXFILENAME];

  int CCDNr;

  /** Projection type (usually SIN). */
  char Projection[MAXMSG];
  /** Right ascension of reference point [deg]. */
  float RefRA;
  /** Declination of reference point [deg]. */
  float RefDec;

  char clobber;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int getpar(struct Parameters* const parameters);


#endif /* ERO_CALEVENTS_H */

