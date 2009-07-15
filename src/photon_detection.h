#ifndef PHOTON_DETECTION_H
#define PHOTON_DETECTION_H 1

#include "sixt.h"
#include "detectors.h"
#include "framestore.h"
#include "photon.h"
#include "eventlistfile.h"
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
  char eventlist_template_filename[FILENAME_LENGTH];
  double t0;
  double timespan;
  
  // Detector specific parameters:
  int detector_type;
  int readout_directions;   // only for DEPFET
  double integration_time;
  double dead_time;
  double clear_time;
  int width; /**< Width of the detector in pixels. */
  double pixelwidth; /**< Width of one detector pixel in [m]. */
  double ccsigma;
  long pha_threshold;
  float energy_threshold;

  float background_rate; /**< Rate of background events. */
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
int getpar(struct Parameters* parameters);


#endif /* PHOTON_DETECTION_H */

