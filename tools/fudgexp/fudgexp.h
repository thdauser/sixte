#ifndef FUDGEXP_H
#define FUDGEXP_H 1

#include "sixt.h"
#include "photon.h"
#include "photonlistfile.h"
#include <wcslib/wcslib.h>

#define TOOLSUB fudgexp_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters{
  /** Filename of the photon list. */
  char PhotonList[MAXFILENAME];

  /** Filename of the exposure map. */
  char ExposureMap[MAXFILENAME];

  /** Nominal exposure time of the map [s]. */
  float ExposureTime;

  char clobber;

  char data_path[MAXFILENAME];
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


#endif /* FUDGEXP_H */
