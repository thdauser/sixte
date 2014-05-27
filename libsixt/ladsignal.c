/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
*/

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

