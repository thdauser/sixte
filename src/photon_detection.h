#ifndef PHOTON_DETECTION_H
#define PHOTON_DETECTION_H 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <limits.h>
#include <assert.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#include "sixt.h"
#include "detectors.h"
#include "photon.h"
#include "eventlist.h"
#include "point.h"


#define TOOLSUB photon_detection_main
#include "headas_main.c"



////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////

struct Parameters {
  char impactlist_filename[FILENAME_LENGTH];
  char attitude_filename[FILENAME_LENGTH];
  char rmf_filename[FILENAME_LENGTH];
  char eventlist_filename[FILENAME_LENGTH];
  double t0;
  double timespan;
  
  // Detector specific parameters:
  int detector_type;
  int readout_directions;   // only for DEPFET
  double integration_time;
  double dead_time;
  double clear_time;
  int width;
  double pixelwidth;
  double ccsigma;
  long pha_threshold;
  float energy_threshold;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int getpar(struct Parameters* parameters);


#endif /* PHOTON_DETECTION_H */

