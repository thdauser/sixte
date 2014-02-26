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

#include "photon.h"


Photon* newPhoton(int* const status)
{
  Photon* ph = (Photon*)malloc(sizeof(Photon));
  CHECK_NULL(ph, *status, "memory allocation for Photon failed");

  // Set initial values.
  ph->ra     = 0.;
  ph->dec    = 0.;
  ph->time   = 0.;
  ph->energy = 0.;
  ph->ph_id  = 0;
  ph->src_id = 0;

  return(ph);
}


void copyPhoton(Photon* const dest, const Photon* const source)
{
  dest->time   = source->time;
  dest->ra     = source->ra;
  dest->dec    = source->dec;
  dest->energy = source->energy;
  dest->ph_id  = source->ph_id;
  dest->src_id = source->src_id;
}

