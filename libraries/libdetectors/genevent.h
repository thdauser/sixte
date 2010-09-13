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

  /** Frame counter. */
  long frame;

  /** Pile-up flag. If the flag is set to 0 the event is not affected
      by pile-up. For any other values the energy (PHA) information
      might be wrong due to energy or pattern pile-up. */
  int pileup;

  /** Split pattern type. Unknown (0), single (1), double (2), triple
      (3), quadruple (4), or invalid (-1). */
  int pat_type;

  /** Pattern ID. Pixel location with the split pattern. */
  int pat_id;
  //   1 2 3
  //   4 5 6
  //   7 8 9

  /** Pattern Alignment. Orientation of the split pattern. */
  int pat_alig;
  //
  // doubles:   2  12  1  21
  //            1      2   
  //        
  // pat_alig:  1  3   5   7
  // -----------------------------------------------------------
  //
  // triples/quaruples:  2   3   12  13  31  21   3   2
  //                     13  12  3   2    2   3  21  31
  //
  // pat_alig:           1   2   3   4    5   6   7   8
  // -----------------------------------------------------------
  // 1-4 : maximum, 2nd, 3rd largest charge.
  
} GenEvent;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////



#endif /* GENEVENT_H */
