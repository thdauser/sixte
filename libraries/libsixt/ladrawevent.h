#ifndef LADRAWEVENT_H 
#define LADRAWEVENT_H 1

#include "sixt.h"


/** Maximum number of photons that are stored as a contribution to a
    single event. If an event originates from more than this
    particular number of photons, the additional ones are not
    stored in the event history. */
#define NLADRAWEVENTPHOTONS (2)


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

  /** Anode ID. */
  long anode;

  /** Event signal. */
  float signal;

  /** Time of event detection [s]. */
  double time;

  /** Identifiers of the contributing photons. */
  long ph_id[NLADRAWEVENTPHOTONS];

  /** Identifiers of the corresponding sources (defined in the SIMPUT
      source catalog). */
  long src_id[NLADRAWEVENTPHOTONS];

} LADRawEvent;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor for LADRawEvent data structure. Initializes pointers
    with NULL and variables with their default values. */
LADRawEvent* getLADRawEvent(int* const status);

/** Destructor for LADRawEvent data structure. */
void freeLADRawEvent(LADRawEvent** const event);


#endif /* LADRAWEVENT_H */
