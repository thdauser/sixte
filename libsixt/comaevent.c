#include "comaevent.h"

CoMaEvent* getCoMaEvent(int* status){
  CoMaEvent* ce=(CoMaEvent*)malloc(sizeof(CoMaEvent));
  if (NULL==ce) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for CoMaEvent!\n",
		   *status);
    return(ce);
  }

  //Initialization:
  ce->time=0.;
  ce->rawx=0;
  ce->rawy=0;
  ce->charge=0.;

  return(ce);
}

void freeCoMaEvent(CoMaEvent** const ce)
{
  if (NULL!=*ce) {
    free(*ce);
    *ce=NULL;
  }
}
