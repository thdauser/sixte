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

  /** Time of event detection. */
  double time;
  
} GenEvent;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////



#endif /* GENEVENT_H */
