#ifndef PHOGEN_H
#define PHOGEN_H 1

#include "sixt.h"

#include "attitudecatalog.h"
#include "geninst.h"
#include "phgen.h"
#include "photonlistfile.h"
#include "simput.h"
#include "vector.h"
#include "simput.h"
#include "sourcecatalog.h"

#define TOOLSUB phogen_main
#include "headas_main.c"


struct Parameters {
  char PhotonList[MAXFILENAME];
  char Mission[MAXMSG];
  char Instrument[MAXMSG];
  char Mode[MAXMSG];
  char XMLFile[MAXFILENAME];
  char Attitude[MAXFILENAME];

  /** [deg] */
  float RA, Dec;

  char Simput[MAXFILENAME];

  double MJDREF;
  double TSTART;
  double Exposure;
  double dt;

  int Seed;
  
  char clobber;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


int phogen_getpar(struct Parameters* parameters);


#endif /* PHOGEN_H */

