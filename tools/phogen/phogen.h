#ifndef PHOGEN_H
#define PHOGEN_H 1

#include "sixt.h"

#include "attitudecatalog.h"
#include "gendet.h"
#include "phgen.h"
#include "photonlistfile.h"
#include "sourceimage.h"
#include "vector.h"
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
  double TIMEZERO;
  double Exposure;

  int Seed;
  
  char clobber;

  char data_path[MAXFILENAME];
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


int phogen_getpar(struct Parameters* parameters);


#endif /* PHOGEN_H */

