#ifndef WFI_SIMULATION_H
#define WFI_SIMULATION_H 1

#include "sixt.h"
#include "detectors.h"
#include "wfidetector.h"
#include "photon.h"
#include "eventlistfile.h"
#include "point.h"
#include "impactlist.h"


#define TOOLSUB wfi_simulation_main
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
  
  // WFI Detector specific parameters:
  int readout_directions;
  double line_readout_time;
  double line_clear_time;

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
static int getpar(struct Parameters* parameters);


#endif /* WFI_SIMULATION_H */

