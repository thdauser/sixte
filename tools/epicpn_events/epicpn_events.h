#ifndef EPICPN_EVENTS_H
#define EPICPN_EVENTS_H 1


#include "sixt.h"
#include "attitude.h"
#include "event.h"
#include "eventfile.h"
#include "wcs.h"

#define TOOLSUB epicpn_events_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters {
  char PatternList[MAXFILENAME];
  char EPICpnEventList[MAXFILENAME];

  char clobber;
};


/** Event entry in an EPIC-pn event file. */
typedef struct {
  /* [s] */
  double time;

  /* [pixel] */
  int rawx, rawy;
  /* [0.05arcsec] */
  int detx, dety;
  /* [0.05arcsec] */
  long x, y;

  /* [adu] */
  int pha;
  /* [eV] */
  int pi;
  
  long flag;
  char pattern;

} EPICpnEvent;


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int getpar(struct Parameters* const parameters);


#endif /* EPICPN_EVENTS_H */

