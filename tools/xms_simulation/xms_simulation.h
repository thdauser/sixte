#ifndef XMS_SIMULATION_H
#define XMS_SIMULATION_H 1

#include "sixt.h"
#include "xmsdetector.h"
#include "eventfile.h"
#include "point.h"
#include "impact.h"
#include "impactlistfile.h"

#define TOOLSUB xms_simulation_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////

struct Parameters {
  char impactlist_filename[MAXFILENAME];
  char rmf_filename_inner[MAXFILENAME], rmf_filename_outer[MAXFILENAME];
  char eventlist_filename[MAXFILENAME];
  char eventlist_template[MAXFILENAME];
  double t0;
  double timespan;
  
  int width_inner, width_outer;
  double pixelwidth_inner, pixelwidth_outer;

  double ccsigma_inner, ccsigma_outer;
  long pha_threshold_inner, pha_threshold_outer;
  float energy_threshold_inner, energy_threshold_outer;

  float background_rate; /**< Rate of background events. */
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
static int getpar(struct Parameters* parameters);


#endif /* XMS_SIMULATION_H */

