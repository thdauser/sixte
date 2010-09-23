#ifndef WFI_PATTERN_RECOMBINATION_H
#define WFI_PATTERN_RECOMBINATION_H 1

#include "sixt.h"
#include "genericdetector.h"
#include "wfidetector.h"

#define TOOLSUB wfi_pattern_recombination_main
#include "headas_main.c"


#define ARRAY_LENGTH 10000


struct Parameters{
  /** Filename of the input WFI event file. */
  char eventlist_filename[MAXMSG];
  /** Filename of the output pattern file. */
  char pattern_filename[MAXMSG];
  /** Template for the new pattern  FITS file. */
  char templatedir[MAXMSG];

  /** Filename of the detector response file containing the EBOUNDS table. */
  char response_filename[MAXMSG];
};


//////////////////////////////////////////////////////////////////


#endif /* WFI_PATTERN_RECOMBINATION_H */
