#ifndef PHOIMG_H
#define PHOIMG_H 1


#include "sixt.h"
#include "attitude.h"
#include "geninst.h"
#include "impactlistfile.h"
#include "photonlistfile.h"
#include "phimg.h"

#define TOOLSUB phoimg_main
#include "headas_main.c"


struct Parameters {
  char PhotonList[MAXFILENAME];
  char ImpactList[MAXFILENAME];
  char Mission[MAXMSG];
  char Instrument[MAXMSG];
  char Mode[MAXMSG];
  char XMLFile[MAXFILENAME];
  char Attitude[MAXFILENAME];

  /** [deg] */
  float RA, Dec;

  double MJDREF;
  double TSTART;
  double Exposure;

  int Seed;
  
  char clobber;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int phoimg_getpar(struct Parameters* parameters);


#endif /* PHOIMG_H */
