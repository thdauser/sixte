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

#ifndef LADSIGNALLIST_H
#define LADSIGNALLIST_H 1

#include "sixt.h"
#include "ladsignal.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Linked, time-ordered list of LADSignals. */
struct structLADSignalListItem {
  LADSignal signal;

  /** Next entry. */
  struct structLADSignalListItem* next;
};
typedef struct structLADSignalListItem LADSignalListItem;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. */
LADSignalListItem* newLADSignalListItem(int* const status);

/** Destructor. */
void freeLADSignalList(LADSignalListItem** const list);


#endif /* LADSIGNALLIST_H */
