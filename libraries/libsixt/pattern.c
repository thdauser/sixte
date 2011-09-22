#include "pattern.h"


Pattern* getPattern(int* const status)
{
  Pattern* pat=(Pattern*)malloc(sizeof(Pattern));
  CHECK_NULL_RET(pat, *status, 
		 "memory allocation for Pattern failed", pat);
  
  // Initalize.
  pat->event=NULL;
  pat->npixels=0;
  pat->type   =0;
  pat->pileup =0;
  int ii;
  for(ii=0; ii<9; ii++) {
    pat->signals[ii]=0.;
  }

  // Call underlying constructors.
  pat->event=getEvent(status);
  CHECK_STATUS_RET(*status, pat);

  return(pat);
}


void freePattern(Pattern** const pattern)
{
  if (NULL!=*pattern) {
    freeEvent(&(*pattern)->event);
    free(*pattern);
    *pattern=NULL;
  }
}
