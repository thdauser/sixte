#ifndef FRAMESTORE_SIMULATION_H
#define FRAMESTORE_SIMULATION_H 1

#include "sixt.h"
#include "detector.h"
#include "framestoredetector.h"
#include "photon.h"
#include "erositaeventfile.h"
#include "point.h"
#include "impactlist.h"


#define TOOLSUB framestore_simulation_main
#include "headas_main.c"



////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////

struct Parameters {
  char impactlist_filename[FILENAME_LENGTH];
  char rmf_filename[MAXMSG];
  char eventlist_filename[MAXMSG];
  char eventlist_template[MAXMSG];
  double t0;
  double timespan;
  
  // Detector specific parameters:
  double integration_time;
  double dead_time;

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


#endif /* FRAMESTORE_SIMULATION_H */

