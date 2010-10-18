#ifndef GENSPLIT_H 
#define GENSPLIT_H 1

#include "sixt.h"
#include "genpixgrid.h"
#include "gendetline.h"
#include "point.h"
#include "genericdetector.h"


// Check only for energy pile-up and ignore pattern pile-up.
//#define ENERGY_PILEUP_ONLY 1

// Check for pattern pile-up in diagonal pixels.
#define DIAGONAL_PATTERN_PILEUP 1


/////////////////////////////////////////////////////////////////
// Constants.
/////////////////////////////////////////////////////////////////


typedef enum {
  GS_NONE        = 0,
  GS_GAUSS       = 1,
  GS_EXPONENTIAL = 2
} GenSplitType;


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic split event generator. */
typedef struct {

  /** Type of the GenSplit model. */
  GenSplitType type;

  /** Split parameter 1. For the Gaussian split model this parameter
      represents the charge cloud size in [m]. For the exponential
      model this parameter represents the denominator in the
      exponential distribution term. */
  double par1;

} GenSplit;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. */
GenSplit* newGenSplit(int* const status);

/** Destructor. */
void destroyGenSplit(GenSplit** const split);

/** Determine split events for a particular photon impact and add the
    fractional charges to the affected pixels. The function return
    value is the number of valid affected pixels inside the detector,
    where charge has been added to. */
int makeGenSplitEvents(const GenSplit* const split,
		       const struct Point2d* const impact,
		       const float charge,
		       const GenPixGrid* const grid,
		       GenDetLine** const line);


#endif /* GENSPLIT_H */
