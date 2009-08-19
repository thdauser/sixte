#ifndef XMS_SIMULATION_H
#define XMS_SIMULATION_H 1

#include "sixt.h"
#include "xmsdetector.h"
#include "eventfile.h"
#include "point.h"
#include "impactlist.h"


#define TOOLSUB xms_simulation_main
#include "headas_main.c"



////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////

struct Parameters {
  char impactlist_filename[FILENAME_LENGTH];
  char rmf_filename_inner[MAXMSG], rmf_filename_outer[MAXMSG];
  char eventlist_filename[MAXMSG];
  char eventlist_template[MAXMSG];
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

