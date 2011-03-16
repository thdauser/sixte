#ifndef GENDETSIM_H
#define GENDETSIM_H 1


#include "sixt.h"
#include "eventlistfile.h"
#include "gendet.h"
#include "impact.h"
#include "impactlistfile.h"

#define TOOLSUB gendetsim_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters {
  char xml_filename[MAXFILENAME];
  char impactlist_filename[MAXFILENAME];
  char eventlist_filename[MAXFILENAME];
  char eventlist_template[MAXFILENAME];

  double t0, exposure;

  int random_seed;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int getpar(struct Parameters* const parameters);


#endif /* GENDETSIM_H */

