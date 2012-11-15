#ifndef COMAIMG_H
#define COMAIMG_H 1


#include "sixt.h"
#include "vector.h"
#include "check_fov.h"
#include "photon.h"
#include "photonlistfile.h"
#include "impact.h"
#include "impactlistfile.h"
#include "telescope.h"
#include "attitudecatalog.h"
#include "codedmask.h"


#define TOOLSUB comaimg_main
#include "headas_main.c"


struct Parameters {
  char Attitude[MAXMSG];   // input: attitude file
  char PhotonList[MAXMSG]; // input: photon list
  char Mask[MAXMSG];       // input: coded mask file
  char ImpactList[MAXMSG]; // output: impact list

  /** [deg] */
  float RA, Dec;

  /** Exposure time. */
  double Exposure;

  // Distance between the coded mask and the detector plane ([m]).
  double MaskDistance;  
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int comaimg_getpar(struct Parameters* parameters);


#endif /* COMAIMG_H */
