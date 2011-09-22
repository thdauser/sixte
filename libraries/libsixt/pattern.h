#ifndef PATTERN_H 
#define PATTERN_H 1

#include "sixt.h"
#include "event.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Event pattern on a pixelized X-ray detector. */
typedef struct {

  Event event;

  /** Split pattern type. Unknown (0), single (1), double (2), triple
      (3), quadruple (4), or invalid (-1). */
  int pat_type;

  /** Pile-up flag. */
  int pileup;

  /** 3x3 array with individual PHA values around the main event. */
  long phas[9];

} Pattern;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////



#endif /* PATTERN_H */
