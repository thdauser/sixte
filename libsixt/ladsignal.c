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


void copyLADSignal(LADSignal* const dest, const LADSignal* const src)
{
  dest->time   = src->time;
  dest->signal = src->signal;
  dest->panel  = src->panel;
  dest->module = src->module;
  dest->element= src->element;
  dest->anode  = src->anode;

  long ii;
  for(ii=0; ii<NLADSIGNALPHOTONS; ii++) {
    dest->ph_id[ii]  = src->ph_id[ii];
    dest->src_id[ii] = src->src_id[ii];
  }
}

