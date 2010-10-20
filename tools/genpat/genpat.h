#ifndef GENPAT_H
#define GENPAT_H 1


#include "sixt.h"
#include "gendet.h"
#include "genevent.h"
#include "geneventfile.h"

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
  char input_eventlist_filename[MAXMSG];
  char output_eventlist_filename[MAXMSG];
};


struct PatternStatistics {
  long nsingles;
  long ndoubles;
  long ntriples;
  long nquadruples;
  long nvalids;
  long ninvalids;
  /** Number of patterns flagged as pile-up. */
  long npileup;
  long npileup_singles;
  /** Number of valid patterns flagged as pile-up. */
  long npileup_valid;
  /** Number of invalid patterns flagged as pile-up. */
  long npileup_invalid;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int getpar(struct Parameters* const parameters);


#endif /* GENPAT_H */

