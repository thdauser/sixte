#ifndef ERO_SPLIT_PHOTONFILE_H
#define ERO_SPLIT_PHOTONFILE_H 1

#include "sixt.h"
#include "photon.h"
#include "photonlistfile.h"

#define TOOLSUB ero_split_photonfile_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters{
  /** Filename of the input photon list. */
  char input_filename[MAXMSG];

  /** Prefix for the output photon list files. */
  char output_prefix[MAXMSG];
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


#endif /* ERO_SPLIT_PHOTONFILE_H */
