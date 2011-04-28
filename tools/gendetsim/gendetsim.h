#ifndef GENDETSIM_H
#define GENDETSIM_H 1


#include "sixt.h"
#include "eventlistfile.h"
#include "gendet.h"
#include "impactlistfile.h"
#include "phdet.h"

#define TOOLSUB gendetsim_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


struct Parameters {
  char ImpactList[MAXFILENAME];
  char EventList[MAXFILENAME];
  char Mission[MAXMSG];
  char Instrument[MAXMSG];
  char Mode[MAXMSG];
  char XMLFile[MAXFILENAME];

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


// Reads the program parameters using PIL
int getpar(struct Parameters* const parameters);


#endif /* GENDETSIM_H */

