#ifndef GENPAT_H
#define GENPAT_H 1


#include "sixt.h"
#include "event.h"
#include "eventlistfile.h"
#include "gendet.h"
#include "genpattern.h"
#include "genpatternfile.h"

#define TOOLSUB genpat_main
#include "headas_main.c"


// Insert only valid events into the output event file with
// the recombined patterns.
#define ONLY_VALID_PATTERNS 1


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters {
  char xml_filename[MAXMSG];
  char eventlist_filename[MAXMSG];
  char patternlist_filename[MAXMSG];
};


struct PatternStatistics {
  /** Number of valid patterns. */
  long nvalids;
  /** Number of valid patterns flagged as pile-up. */
  long npvalids;

  /** Number of invalid patterns. */
  long ninvalids;
  /** Number of invalid patterns flagged as pile-up. */
  long npinvalids;

  /** Number of patterns with a particular grade. */
  long ngrade[256];
  /** NUmber of patterns with a particular grade flagged as pile-up. */
  long npgrade[256];
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int getpar(struct Parameters* const parameters);


#endif /* GENPAT_H */

