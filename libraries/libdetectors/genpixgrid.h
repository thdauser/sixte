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
void destroyGenPixGrid(GenPixGrid** const grid);

/** Return the index of the detector line affected by the specified
    y-value. If the y-value is outside the line region, the return
    value is -1. */
inline int getGenDetAffectedLine(const GenPixGrid* const grid, const double y);

/** Return the index of the detector column affected by the specified
    x-value. If the x-value is outside the line region, the return
    value is -1. */
inline int getGenDetAffectedColumn(const GenPixGrid* const grid, const double x);


#endif /* GENPIXGRID_H */
