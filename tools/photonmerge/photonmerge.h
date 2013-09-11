#ifndef PHOTONMERGE_H
#define PHOTONMERGE_H 1

#include "sixt.h"
#include "photon.h"
#include "photonfile.h"

#define TOOLSUB photonmerge_main
#include "headas_main.c"

#define MAX_N_INPUTFILES 1000

////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters{
  /** Number of input files. */
  int NInputFiles;

  /** Filename prefix of the input photon list files. */
  char InputPrefix[MAXFILENAME];

  /** Filename of the output photon list file. */
  char OutputFile[MAXFILENAME];

  char clobber;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


#endif /* PHOTONMERGE_H */
