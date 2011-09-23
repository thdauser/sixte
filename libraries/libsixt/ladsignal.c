#include "ladsignal.h"


LADSignal* getLADSignal(int* const status)
{
  LADSignal* ev = (LADSignal*)malloc(sizeof(LADSignal));
  CHECK_NULL_RET(ev, *status, 
		 "memory allocation for LADSignal failed", ev);
  
  // Initalize.
  ev->panel  =0;
  ev->module =0;
  ev->element=0;
  ev->anode  =0;
  ev->signal =0.;
  ev->time   =0.;

  long ii;
  for(ii=0; ii<NLADSIGNALPHOTONS; ii++) {
    ev->ph_id[ii]  = 0;
    ev->src_id[ii] = 0;
  }

  return(ev);
}


void freeLADSignal(LADSignal** const event)
{
  if (NULL!=*event) {
    free(*event);
    *event=NULL;
  }
}

