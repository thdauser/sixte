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

#ifndef LADIMPACT_H
#define LADIMPACT_H 1

#include "sixt.h"
#include "point.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic event on a pixelized X-ray detector. */
typedef struct {

  /** Panel ID. */
  long panel;

  /** Module ID. */
  long module;

  /** Element ID. */
  long element;

  /** Impact position of the photon on the element [m]. */
  struct Point2d position;

  /** Photon energy [keV]. */
  float energy;

  /** Arrival time of the photon [s]. */
  double time;

  /** Identifier of the photon. */
  long ph_id;

  /** Identifier of the corresponding source (defined in the SIMPUT
      source catalog). */
  long src_id;

} LADImpact;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor for LADImpact data structure. Initializes pointers with
    NULL and variables with their default values. */
LADImpact* getLADImpact(int* const status);

/** Destructor for LADImpact data structure. */
void freeLADImpact(LADImpact** const impact);


#endif /* LADIMPACT_H */
