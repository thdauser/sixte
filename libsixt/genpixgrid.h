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
  /** Pixel width [m]. The pixel width includes (2 times) the pixel border. */
  float xdelt, ydelt;
  /** Pixel border [m]. In this area along the corner the pixel is
      insensitive to incident photons. If a photon hits this border
      area, the return value of getGenDetAffectedLine() or
      getGenDetAffectedColumn() is -1. */
  float xborder, yborder;

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
int getGenDetAffectedLine(const GenPixGrid* const grid, const double y);

/** Return the index of the detector column affected by the specified
    x-value. If the x-value is outside the line region, the return
    value is -1. */
int getGenDetAffectedColumn(const GenPixGrid* const grid, const double x);


#endif /* GENPIXGRID_H */
