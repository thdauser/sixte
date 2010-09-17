#ifndef GENPAT_H
#define GENPAT_H 1


#include "sixt.h"
#include "gendet.h"
#include "genevent.h"
#include "geneventfile.h"

#define TOOLSUB genpat_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters {
  char xml_filename[MAXMSG];
  char input_eventlist_filename[MAXMSG];
  char output_eventlist_filename[MAXMSG];
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int getpar(struct Parameters* const parameters);


#endif /* GENPAT_H */

