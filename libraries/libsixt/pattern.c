#include "pattern.h"


Pattern* getPattern(int* const status)
{
  Pattern* pat=(Pattern*)malloc(sizeof(Pattern));
  CHECK_NULL_RET(pat, *status, 
		 "memory allocation for Pattern failed", pat);
  
  // Initalize.
  pat->rawx=0;
  pat->rawy=0;
  pat->pha =0;
  pat->signal =0.;
  pat->time   =0.;
  pat->frame  =0;
  pat->npixels=0;
  pat->type   =0;
  pat->pileup =0;
  pat->ra     =0.;
  pat->dec    =0.;

  int ii;
  for(ii=0; ii<9; ii++) {
    pat->signals[ii]=0.;
  }
  for(ii=0; ii<NPATTERNPHOTONS; ii++) {
    pat->ph_id[ii]  = 0;
    pat->src_id[ii] = 0;
  }

  return(pat);
}


void freePattern(Pattern** const pattern)
{
  if (NULL!=*pattern) {
    free(*pattern);
    *pattern=NULL;
  }
}
