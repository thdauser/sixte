#ifndef PHOIMG_H
#define PHOIMG_H 1


#include "sixt.h"
#include "attitudecatalog.h"
#include "gendet.h"
#include "impactlistfile.h"
#include "photonlistfile.h"
#include "phimg.h"

#define TOOLSUB phoimg_main
#include "headas_main.c"


struct Parameters {
  char attitude_filename[MAXMSG];   // input: attitude file
  char photonlist_filename[MAXMSG]; // input: photon list
  char xml_filename[MAXMSG];        // input: detector XML description
  char impactlist_filename[MAXMSG]; // output: impact list
  char impactlist_template[MAXMSG];

  double t0;       // starting time of the simulation
  double exposure; // time span of the simulation

  int random_seed;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int phoimg_getpar(struct Parameters* parameters);


#endif /* PHOIMG_H */
