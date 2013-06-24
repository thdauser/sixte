#ifndef ERO_RAWEVENTS_H
#define ERO_RAWEVENTS_H 1


#include "sixt.h"
#include "event.h"
#include "eventlistfile.h"

#define TOOLSUB ero_rawevents_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters {
  char EventList[MAXFILENAME];
  char eroEventList[MAXFILENAME];

  char clobber;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int getpar(struct Parameters* const parameters);


#endif /* ERO_RAWEVENTS_H */

