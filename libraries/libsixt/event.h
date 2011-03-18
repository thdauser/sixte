#ifndef EVENT_H 
#define EVENT_H 1

#include "sixt.h"


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

} Event;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////



#endif /* EVENT_H */
