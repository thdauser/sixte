#ifndef PHOIMG_H
#define PHOIMG_H 1


#include "sixt.h"
#include "attitudecatalog.h"
#include "gendet.h"
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
  double TIMEZERO;
  double Exposure;

  int Seed;
  
  char clobber;

  char fits_templates[MAXFILENAME];
  char xml_path[MAXFILENAME];
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int phoimg_getpar(struct Parameters* parameters);


#endif /* PHOIMG_H */
