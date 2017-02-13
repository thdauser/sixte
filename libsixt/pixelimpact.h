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

#ifndef PIXIMPACT_H
#define PIXIMPACT_H 1

#include "sixt.h"
#include "point.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Impact of a photon on the detector plane. */
typedef struct {
  
  /** ID of pixel. */
  long pixID;
  
  /** Arrival time of the photon on the detector [s]. */
  double time;
  
  /** Photon energy [keV]. */
  float energy;
  
  /** Impact position of the photon on the detector [m]. */
  struct Point2d detposition;
  
  /** Impact position of the photon on the pixel [m]. */
  struct Point2d pixposition;

  /** Unique photon identifier. */
  long ph_id;

  /** Unique source identifier for the originating X-ray source. */
  long src_id;

  /** time in samples after this impact **/
  long grade1;

  /** time in samples before this impact **/
  long grade2;

  /** total energy of the final impact potentially containing energy from previous impact (pileup) */
  double totalenergy;

  /** number of parent events */
  int nb_pileup;

  /** weight index (for ctk) */
  int weight_index;

} PixImpact;


#endif /* PIXIMPACT_H */
