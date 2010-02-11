#ifndef PNCCD_SIMULATION_H
#define PNCCD_SIMULATION_H

#include "sixt.h"
#include "pnccddetector.h"
#include "eventfile.h"
#include "point.h"
#include "impact.h"
#include "impactlistfile.h"

#define TOOLSUB pnccd_simulation_main
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
  
  double ccsigma;
  long pha_threshold;
  float energy_threshold;

  float background_rate; /**< Rate of background events. */

  // HTRS Detector specific parameters:
  double pixelwidth;
  double dead_time;
};
