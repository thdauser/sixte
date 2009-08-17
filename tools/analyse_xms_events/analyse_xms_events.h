#ifndef ANALYSE_XMS_EVENTS_H
#define ANALYSE_XMS_EVENTS_H 1

#include "sixt.h"
#include "xmseventfile.h"


#define TOOLSUB analyse_xms_events_main
#include "headas_main.c"

// Number of photons stored in the event list for each individual pixel.
#define XMS_EVENT_LIST_ENTRIES 5


struct Parameters{
  /** Filename of the XMS event file. */
  char eventlist_filename[MAXMSG];

  /** Characteristic time unit of the TES microcalorimeter. */
  double time_unit;
  /** Time units before and after and event that may not be affected by further
   * impacting photons in the same pixel. Otherwise the event energy determination
   * will be degraded. */
  int units_before_pulse, units_after_pulse;
};


//////////////////////////////////////////////////////////////////


int analyse_xms_events_getpar(struct Parameters*);


#endif /* ANALYSE_XMS_EVENTS_H */
