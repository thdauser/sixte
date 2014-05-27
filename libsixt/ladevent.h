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

#ifndef LADEVENT_H 
#define LADEVENT_H 1

#include "sixt.h"


/** Maximum number of photons that are stored as a contribution to a
    single event. If an event originates from more than this
    particular number of photons, the additional ones are not
    stored in the event history. */
#define NLADEVENTPHOTONS (2)


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

  /** Anode ID of center anode. */
  long anode;

  /** Event signal [keV]. */
  float signal;

  /** Time of event detection [s]. */
  double time;

  /** Identifiers of the contributing photons. */
  long ph_id[NLADEVENTPHOTONS];

  /** Identifiers of the corresponding sources (defined in the SIMPUT
      source catalog). */
  long src_id[NLADEVENTPHOTONS];

} LADEvent;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor for LADEvent data structure. Initializes pointers
    with NULL and variables with their default values. */
LADEvent* getLADEvent(int* const status);

/** Destructor for LADEvent data structure. */
void freeLADEvent(LADEvent** const event);

/** Copy LADEvent data structure. */
void copyLADEvent(LADEvent* const dest, const LADEvent* const src);


#endif /* LADEVENT_H */
