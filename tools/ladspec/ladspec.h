#ifndef LADSPEC_H
#define LADSPEC_H 1

#include "sixt.h"
#include "lad.h"
#include "ladeventlistfile.h"
#include "ladevent.h"

#define TOOLSUB ladspec_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char EventList[MAXFILENAME];
  char Spectrum[MAXFILENAME];
  char XMLFile[MAXFILENAME];

  char clobber;
};


int ladspec_getpar(struct Parameters *par);


#endif /* LADSPEC_H */

