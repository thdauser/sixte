#ifndef GENEVENT_H 
#define GENEVENT_H 1

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

  /** Pile-up flag. If the flag is set to 0 the event is not affected
      by pile-up. For any other values the energy (PHA) information
      might be wrong due to energy or pattern pile-up. */
  int pileup;
  
} GenEvent;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////



#endif /* GENEVENT_H */
