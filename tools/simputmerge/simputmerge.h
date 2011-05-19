#ifndef SIMPUTMERGE_H
#define SIMPUTMERGE_H 1

#include "sixt.h"
#include "simput.h"

#define TOOLSUB simputmerge_main
#include "headas_main.c"


struct Parameters {
  char Infile1[MAXFILENAME];
  char Infile2[MAXFILENAME];
  char Outfile[MAXFILENAME];
  char FetchExtensions;
  
  char clobber;
};


int simputmerge_getpar(struct Parameters* const par);


#endif /* SIMPUTMERGE_H */

