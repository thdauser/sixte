#ifndef EVENT_H 
#define EVENT_H 1

#include "sixt.h"


/** Maximum number of photons that are stored as a contribution to a
    single event. If an event originates from more than this
    particular number of photons, the additional ones are not stored
    in the event history. */
#define NEVENTPHOTONS (2)


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Event on a pixelized X-ray detector. */
typedef struct {

  /** Raw detector coordinates. Indices start at 0. */
  int rawx, rawy;

  /** Detected PI channel [adu]. */
  long pi;

  /** Signal in [keV]. */
  float signal;

  /** Time of detection [s]. */
  double time;

  /** Frame counter. */
  long frame;
  
  /** Back-projected right ascension to the sky [rad]. */
  double ra;

  /** Back-projected declination to the sky [rad]. */
  double dec;

  /** Identifiers of the contributing photons. */
  long ph_id[NEVENTPHOTONS];

  /** Identifiers of the corresponding sources (defined in the SIMPUT
      source catalog). */
  long src_id[NEVENTPHOTONS];

  /** Number of pixels involved in the split pattern. This number is
      one for a single-pixel event. */
  long npixels;

  /** Split pattern type: single (0), double (1-4), triple (5-8),
      quadruple (9-12), or invalid (-1). */
  int type;

  /** Pile-up flag. */
  int pileup;

  /** 3x3 array with individual signal values [keV] around the main
      event. */
  float signals[9];

  /** 3x3 array with individual energy channels [adu] around the main
      event. */
  long pis[9];

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
