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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef IMPACT_H
#define IMPACT_H 1

#include "sixt.h"
#include "point.h"

/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Impact of a photon on the detector plane. */
typedef struct {
  /** Arrival time of the photon on the detector [s]. */
  double time;

  /** Photon energy [keV]. */
  float energy;

  /** Impact position of the photon on the detector [m]. */
  struct Point2d position;

  /** Unique photon identifier. */
  long ph_id;

  /** Unique source identifier for the originating X-ray source. */
  long src_id;

} Impact;


Impact* newImpact(int* const status);

void copyImpact(Impact* const dest, const Impact* const source);

#endif /* IMPACT_H */
