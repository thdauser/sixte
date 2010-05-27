#ifndef ERO_MERGE_PHOTONFILES_H
#define ERO_MERGE_PHOTONFILES_H 1

#include "sixt.h"
#include "photon.h"
#include "photonlistfile.h"

#define TOOLSUB ero_merge_photonfiles_main
#include "headas_main.c"

#define MAX_N_INPUTFILES 100

////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters{
  /** Number of input files. */
  int n_inputfiles;

  /** Filename prefix of the input event files. */
  char input_prefix[MAXMSG];

  /** Filename of the output event file. */
  char output_filename[MAXMSG];

  char photonlist_template[MAXMSG];
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


#endif /* ERO_MERGE_PHOTONFILES_H */
