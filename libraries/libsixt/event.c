#include "event.h"


Event* getEvent(int* const status)
{
  Event* ev = (Event*)malloc(sizeof(Event));
  CHECK_NULL_RET(ev, *status, 
		 "memory allocation for Event failed", ev);
  
  // Initalize.
  ev->rawx=0;
  ev->rawy=0;
  ev->pha=0;
  ev->charge=0.;
  ev->time=0.;
  ev->frame=0;

  int ii;
  for(ii=0; ii<NEVENTPHOTONS; ii++) {
    ev->ph_id[ii]  = 0;
    ev->src_id[ii] = 0;
  }

  return(ev);
}


void freeEvent(Event** const event)
{
  if (NULL!=*event) {
    free(*event);
    *event=NULL;
  }
}
