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

#include "event.h"


Event* getEvent(int* const status)
{
  Event* ev=(Event*)malloc(sizeof(Event));
  CHECK_NULL_RET(ev, *status, 
		 "memory allocation for Event failed", ev);
  
  // Initalize.
  ev->rawx=0;
  ev->rawy=0;
  ev->pha  =0;
  ev->signal =0.;
  ev->time   =0.;
  ev->frame  =0;
  ev->npixels=0;
  ev->type   =0;
  ev->pileup =0;
  ev->ra     =0.;
  ev->dec    =0.;

  int ii;
  for(ii=0; ii<9; ii++) {
    ev->signals[ii]=0.;
    ev->phas[ii]    =0;
  }
  for(ii=0; ii<NEVENTPHOTONS; ii++) {
    ev->ph_id[ii] =0;
    ev->src_id[ii]=0;
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
