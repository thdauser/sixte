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

#ifndef LINKEDIMPLIST_H
#define LINKEDIMPLIST_H 1

#include "sixt.h"
#include "impact.h"

/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Linked, time-ordered photon list. */
struct structLinkedImpListElement {

  Impact impact;

  /** Next entry. */
  struct structLinkedImpListElement* next;
};
typedef struct structLinkedImpListElement LinkedImpListElement;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////

/** Constructor. */
LinkedImpListElement* newLinkedImpListElement(int* const status);

/** Destructor. */
void freeLinkedImpList(LinkedImpListElement** const list);

/** Merge 2 time-ordered linked photon lists. */
LinkedImpListElement* mergeLinkedImpLists(LinkedImpListElement* list1,
					  LinkedImpListElement* list2);


#endif /* LINKEDIMPLIST_H */
