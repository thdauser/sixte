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
	// Names of the used files
  char impactlist_filename[FILENAME_LENGTH];
  char rmf_filename[MAXMSG];
  char eventlist_filename[MAXMSG];
  char eventlist_template[MAXMSG];
	// Times needed for the simulation
  double t0;
  double timespan;
  
	// Detector parameters
  double ccsigma; /* Sigma of the charge cloud. Has to be maybe changed due to the king profile for XMM */
  long pha_threshold;
  float energy_threshold;

	// Source/Detector parameter
  float background_rate; /**< Rate of background events. */

	// Detector
	int xwidth;
	int ywidth;
  double pixelwidth; /**< Width of one detector pixel in [m]. */
  double dead_time;
	
	// Readout
	char readout_mode[MAXMSG];
	int readout_directions;
	double readout_time;
	double clear_time;

};

// Reads the program parameters using PIL
static int getpar(struct Parameters* parameters);

#endif /* PNCCD_SIMULATION_H */
