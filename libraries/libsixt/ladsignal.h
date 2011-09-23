#ifndef LADSIGNAL_H 
#define LADSIGNAL_H 1

#include "sixt.h"


/** Maximum number of photons that are stored as a contribution to a
    single event. If an event originates from more than this
    particular number of photons, the additional ones are not
    stored in the event history. */
#define NLADSIGNALPHOTONS (2)


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

  /** Event signal [keV]. */
  float signal;

  /** Time of event detection [s]. */
  double time;

  /** Identifiers of the contributing photons. */
  long ph_id[NLADSIGNALPHOTONS];

  /** Identifiers of the corresponding sources (defined in the SIMPUT
      source catalog). */
  long src_id[NLADSIGNALPHOTONS];

} LADSignal;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor for LADSignal data structure. Initializes pointers
    with NULL and variables with their default values. */
LADSignal* getLADSignal(int* const status);

/** Destructor for LADSignal data structure. */
void freeLADSignal(LADSignal** const event);


#endif /* LADSIGNAL_H */
