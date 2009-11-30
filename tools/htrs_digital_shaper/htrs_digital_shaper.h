#ifndef HTRS_DIGITAL_SHAPER_H
#define HTRS_DIGITAL_SHAPER_H 1

#include "sixt.h"
#include "htrseventfile.h"
#include "htrsevent.h"


#define TOOLSUB htrs_digital_shaper_main
#include "headas_main.c"


struct Parameters{
  /** Filename of the input HTRS event file containing all events. */
  char input_eventlist_filename[MAXMSG];
  /** Filename of the output HTRS event file containing only events
      that are properly measured by the digital shaper. */
  char output_eventlist_filename[MAXMSG];
  /** Filename of the template for a new HTRS event file. */
  char eventlist_template[MAXMSG];

  /** Shaping frequency (Hz). */
  double frequency;
  /** Number of sampling points. */
  int nsamplings;
};


//////////////////////////////////////////////////////////////////


int htrs_digital_shaper_getpar(struct Parameters*);


#endif /* HTRS_DIGITAL_SHAPER_H */
