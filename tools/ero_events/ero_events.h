#ifndef ERO_EVENTS_H
#define ERO_EVENTS_H 1


#include "sixt.h"
#include "event.h"
#include "eventlistfile.h"

#define TOOLSUB ero_events_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters {
  char genEventList[MAXFILENAME];
  char eroEventList[MAXFILENAME];

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


#endif /* ERO_EVENTS_H */

