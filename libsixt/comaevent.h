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

#ifndef COMAEVENT_H
#define COMAEVENT_H 1

#include "sixt.h"

/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////

/** Coded Mask specific event. */
typedef struct {
  double time;
  /** Pixel charge represented by the photon energy [keV]. */
  double charge; 
  /** Pixel coordinates starting at 0. */
  int rawx, rawy;
} CoMaEvent;

/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////

CoMaEvent* getCoMaEvent(int* status);

void freeCoMaEvent(CoMaEvent** const ce);


#endif /* COMAEVENT_H */

