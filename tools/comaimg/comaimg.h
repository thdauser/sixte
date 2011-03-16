#ifndef COMAIMG_H
#define COMAIMG_H 1


#include "sixt.h"
#include "vector.h"
#include "check_fov.h"
#include "photon.h"
#include "photonlistfile.h"
#include "impact.h"
#include "impactlistfile.h"
#include "telescope.h"
#include "attitudecatalog.h"
#include "codedmask.h"


#define TOOLSUB comaimg_main
#include "headas_main.c"


struct Parameters {
  char attitude_filename[MAXMSG];   // input: attitude file
  char photonlist_filename[MAXMSG]; // input: photon list
  char mask_filename[MAXMSG];       // input: coded mask file
  char impactlist_filename[MAXMSG]; // output: impact list
  char impactlist_template[MAXMSG]; // template for the impact list FITS file

  // Distance between the coded mask and the detector plane ([m]).
  double mask_distance;  
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int comaimg_getpar(struct Parameters* parameters);


#endif /* COMAIMG_H */
