#ifndef GENPIXGRID_H 
#define GENPIXGRID_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic square pixel grid. */
typedef struct {

  /** Detector dimensions. Width and height [pixels]. */
  int xwidth, ywidth;

  /** Reference pixel. */
  float xrpix, yrpix;
  /** Reference value [m]. */
  float xrval, yrval;
  /** Pixel width [m]. */
  float xdelt, ydelt;

} GenPixGrid;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. */
GenPixGrid* newGenPixGrid(int* const status);

/** Destructor. */
void destroyGenPixGrid(GenPixGrid** grid);


#endif /* GENPIXGRID_H */
