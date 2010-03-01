#ifndef ERO_ONBOARD_PROC_H
#define ERO_ONBOARD_PROC_H 1

#include "sixt.h"
#include "genericdetector.h"
#include "framestoredetector.h"

#define TOOLSUB ero_onboard_proc_main
#include "headas_main.c"


#define ARRAY_LENGTH 1000


struct Parameters{
  /** Filename of the input eROSITA event file. */
  char eventlist_filename[MAXMSG];

  /** Filename of the detector response file containing the EBOUNDS table. */
  char response_filename[MAXMSG];

  long pha_threshold;
  float energy_threshold;

};


//////////////////////////////////////////////////////////////////


#endif /* ERO_ONBOARD_PROC_H */
