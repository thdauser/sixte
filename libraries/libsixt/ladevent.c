#include "ladevent.h"


LADEvent* getLADEvent(int* const status)
{
  LADEvent* ev = (LADEvent*)malloc(sizeof(LADEvent));
  CHECK_NULL_RET(ev, *status, 
		 "memory allocation for LADEvent failed", ev);
  
  // Initalize.
  ev->panel  =0;
  ev->module =0;
  ev->element=0;
  ev->anode  =0;
  ev->signal =0.;
  ev->time   =0.;

  long ii;
  for(ii=0; ii<NLADEVENTPHOTONS; ii++) {
    ev->ph_id[ii]  = 0;
    ev->src_id[ii] = 0;
  }

  return(ev);
}


void freeLADEvent(LADEvent** const event)
{
  if (NULL!=*event) {
    free(*event);
    *event=NULL;
  }
}

