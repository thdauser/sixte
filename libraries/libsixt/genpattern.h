#ifndef GENPATTERN_H 
#define GENPATTERN_H 1

#include "sixt.h"
#include "event.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic event pattern on a pixelized X-ray detector. */
typedef struct {

  Event event;

  /** Split pattern type. Unknown (0), single (1), double (2), triple
      (3), quadruple (4), or invalid (-1). */
  int pat_type;

  /** 3x3 array with individual PHA values around the main event. */
  long phas[9];

} GenPattern;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////



#endif /* GENPATTERN_H */
