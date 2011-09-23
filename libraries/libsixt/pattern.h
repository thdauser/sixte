#ifndef PATTERN_H 
#define PATTERN_H 1

#include "sixt.h"


/** Maximum number of photons that are stored as a contribution to a
    single pattern. If a pattern originates from more than this
    particular number of photons, the additional ones are not stored
    in the event history. */
#define NPATTERNPHOTONS (2)


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Event pattern on a pixelized X-ray detector. */
typedef struct {

  /** Raw detector coordinates. Indices start at 0. */
  int rawx, rawy;

  /** Detected PHA channel. */
  long pha;

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
  long ph_id[NPATTERNPHOTONS];

  /** Identifiers of the corresponding sources (defined in the SIMPUT
      source catalog). */
  long src_id[NPATTERNPHOTONS];

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
