#ifndef SIMPUTPSD_H
#define SIMPUTPSD_H 1

#include "sixt.h"
#include "simput.h"

#define TOOLSUB simputpsd_main
#include "headas_main.c"


struct Parameters {
  /** File name of the SIMPUT file, where the PSD should be attached
      to. */
  char Simput[MAXFILENAME];
  
  /** File name of the input ASCII PSD. */
  char PSDFile[MAXFILENAME];
};


int simputpsd_getpar(struct Parameters* const par);


#endif /* SIMPUTPSD_H */

