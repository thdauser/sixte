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

  char clobber;

  char fits_templates[MAXFILENAME];
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int getpar(struct Parameters* const parameters);


#endif /* ERO_EVENTS_H */

