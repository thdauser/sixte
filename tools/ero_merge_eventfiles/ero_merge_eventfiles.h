#ifndef ERO_MERGE_EVENTFILES_H
#define ERO_MERGE_EVENTFILES_H 1

#include "sixt.h"
#include "erositaevent.h"
#include "erositaeventfile.h"

#define TOOLSUB ero_merge_eventfiles_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters{
  /** Filename of the input event files. */
  char input_prefix[MAXMSG];

  /** Filename of the output event file. */
  char output_filename[MAXMSG];

  char eventlist_template[MAXMSG];
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


#endif /* ERO_MERGE_EVENTFILES_H */
