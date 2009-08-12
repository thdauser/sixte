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
  char rmf_filename[MAXMSG];
  char eventlist_filename[MAXMSG];
  char eventlist_template[MAXMSG];
  double t0;
  double timespan;
  
  int width;
  double pixelwidth;

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


#endif /* XMS_SIMULATION_H */

