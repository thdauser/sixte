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
