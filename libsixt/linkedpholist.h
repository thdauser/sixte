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

#ifndef LINKEDPHOLIST_H
#define LINKEDPHOLIST_H 1

#include "sixt.h"
#include "photon.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Linked, time-ordered photon list. */
struct structLinkedPhoListElement {

  Photon photon;

  /** Next entry. */
  struct structLinkedPhoListElement* next;
};
typedef struct structLinkedPhoListElement LinkedPhoListElement;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. */
LinkedPhoListElement* newLinkedPhoListElement(int* const status);

/** Destructor. */
void freeLinkedPhoList(LinkedPhoListElement** const list);

/** Merge 2 time-ordered linked photon lists. */
LinkedPhoListElement* mergeLinkedPhoLists(LinkedPhoListElement* list1, 
					  LinkedPhoListElement* list2);


#endif /* LINKEDPHOLIST_H */
