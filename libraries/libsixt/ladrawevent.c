#include "ladrawevent.h"


LADRawEvent* getLADRawEvent(int* const status)
{
  LADRawEvent* ev = (LADRawEvent*)malloc(sizeof(LADRawEvent));
  CHECK_NULL_RET(ev, *status, 
		 "memory allocation for LADRawEvent failed", ev);
  
  // Initalize.
  ev->panel  =0;
  ev->module =0;
  ev->element=0;
  ev->anode  =0;
  ev->signal =0.;
  ev->time   =0.;

  long ii;
  for(ii=0; ii<NLADRAWEVENTPHOTONS; ii++) {
    ev->ph_id[ii]  = 0;
    ev->src_id[ii] = 0;
  }

  return(ev);
}


void freeLADRawEvent(LADRawEvent** const event)
{
  if (NULL!=*event) {
    free(*event);
    *event=NULL;
  }
}

