#ifndef PATTERN_H 
#define PATTERN_H 1

#include "sixt.h"
#include "event.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Event pattern on a pixelized X-ray detector. */
typedef struct {

  Event* event;

  /** Number of pixels involved in the pattern. */
  long npixels;

  /** Split pattern type. Unknown (0), single (1), double (2), triple
      (3), quadruple (4), or invalid (-1). */
  int type;

  /** Pile-up flag. */
  int pileup;

  /** 3x3 array with individual signal values [keV] around the main
      event. */
  float signals[9];

} Pattern;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor for Pattern data structure. Initializes pointers with
    NULL and variables with their default values. */
Pattern* getPattern(int* const status);

/** Destructor for Pattern data structure. */
void freePattern(Pattern** const pattern);


#endif /* PATTERN_H */
