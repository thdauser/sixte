#ifndef ERO_CALEVENTS_H
#define ERO_CALEVENTS_H 1


#include "sixt.h"
#include "attitude.h"
#include "event.h"
#include "eventfile.h"
#include "gti.h"
#include "wcs.h"

#define TOOLSUB ero_calevents_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters {
  char PatternList[MAXFILENAME];
  char eroEventList[MAXFILENAME];
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


/** Event entry in an eROSITA calibrated event file. */
typedef struct {
  double time; /* [s] */
  long frame;

  long pha; /* [adu] */
  float energy; /* [eV] */

  int rawx, rawy;
  int subx, suby;

  double ra, dec; /* [deg] */
  long x, y;

  long flag;

  unsigned int pat_typ;
  unsigned char pat_inf;

  float ev_weight;

  unsigned char ccdnr;

} eroCalEvent;


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int getpar(struct Parameters* const parameters);


#endif /* ERO_CALEVENTS_H */

