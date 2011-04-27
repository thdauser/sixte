#ifndef EVENT_H 
#define EVENT_H 1

#include "sixt.h"


/** Maximum number of photons that are stored as a contribution to a
    single event. If an event originates from more than this
    particular number of photons, the additional ones are not
    stored in the event history. */
#define NEVENTPHOTONS (2)


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic event on a pixelized X-ray detector. */
typedef struct {
  
  /** Raw detector coordinates. Indices start at 0. */
  int rawx, rawy;

  /** Detected PHA channel. */
  long pha;

  /** Pixel charge in [keV]. */
  float charge;

  /** Time of event detection. */
  double time;

  /** Frame counter. */
  long frame;

  /** Identifiers of the contributing photons. */
  long ph_id[NEVENTPHOTONS];

  /** Identifiers of the corresponding sources (defined in the SIMPUT
      source catalog). */
  long src_id[NEVENTPHOTONS];

} Event;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor for Event data structure. Initializes pointers with
    NULL and variables with their default values. */
Event* getEvent(int* const status);

/** Destructor for Event data structure. */
void freeEvent(Event** const event);


#endif /* EVENT_H */
