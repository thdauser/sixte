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
