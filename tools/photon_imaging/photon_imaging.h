#ifndef PHOTON_IMAGING_H
#define PHOTON_IMAGING_H 1


#include "sixt.h"
#include "vector.h"
#include "check_fov.h"
#include "sixt_random.h"
#include "psf.h"
#include "photon.h"
#include "photonlistfile.h"
#include "impact.h"
#include "impactlistfile.h"
#include "telescope.h"
#include "attitudecatalog.h"
#include "vignetting.h"


#define TOOLSUB photon_imaging_main
#include "headas_main.c"


struct Parameters {
  char attitude_filename[MAXMSG];   // input: attitude file
  char photonlist_filename[MAXMSG]; // input: photon list
  char psf_filename[MAXMSG];        // input: PSF file
  char vignetting_filename[MAXMSG]; // input: vignetting file
  char impactlist_filename[MAXMSG]; // output: impact list
  char impactlist_template[MAXMSG];

  /** Telescope number */
  int telescope;

  double t0;        // starting time of the simulation
  double timespan;  // time span of the simulation
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int photon_imaging_getpar(struct Parameters* parameters, struct Telescope *);


#endif /* PHOTON_IMAGING_H */
