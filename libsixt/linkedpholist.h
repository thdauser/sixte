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
