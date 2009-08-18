#ifndef PATTERN_RECOMBINATION_H
#define PATTERN_RECOMBINATION_H 1

#include "sixt.h"
#include "wfidetector.h"

#define TOOLSUB pattern_recombination_main
#include "headas_main.c"


struct Parameters{
  /** Filename of the XMS event file. */
  char eventlist_filename[MAXMSG];
};


//////////////////////////////////////////////////////////////////


int pattern_recombination_getpar(struct Parameters*);


#endif /* PATTERN_RECOMBINATION_H */
