#ifndef PHOTON_IMAGING_H
#define PHOTON_IMAGING_H 1


#include "sixt.h"
#include "vector.h"
#include "check_fov.h"
#include "random_sixt.h"
#include "psf.h"
#include "photon.h"
#include "telescope.h"
#include "attitudecatalog.h"
#include "vignetting.h"


#define TOOLSUB photon_imaging_main
#include "headas_main.c"


struct Parameters {
  char attitude_filename[FILENAME_LENGTH];
  char photonlist_filename[FILENAME_LENGTH]; // input: photon list
  char psf_filename[FILENAME_LENGTH];        // input: PSF input file
  char vignetting_filename[FILENAME_LENGTH];
  char impactlist_filename[FILENAME_LENGTH]; // output: impact list

  double t0;        // starting time of the simulation
  double timespan;  // time span of the simulation
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int photon_imaging_getpar(struct Parameters* parameters, struct Telescope *);


#endif /* PHOTON_IMAGING_H */
