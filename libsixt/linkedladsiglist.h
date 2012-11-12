#ifndef LINKEDLADSIGLIST_H
#define LINKEDLADSIGLIST_H 1

#include "sixt.h"
#include "ladsignal.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Linked, time-ordered list of LADSignals. */
struct structLinkedLADSigListElement {

  LADSignal signal;

  /** Next entry. */
  struct structLinkedLADSigListElement* next;
};
typedef struct structLinkedLADSigListElement LinkedLADSigListElement;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. */
LinkedLADSigListElement* newLinkedLADSigListElement(int* const status);

/** Destructor. */
void freeLinkedLADSigList(LinkedLADSigListElement** const list);


#endif /* LINKEDLADSIGLIST_H */
