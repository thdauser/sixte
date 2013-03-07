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

  /** Pixel width [m]. The pixel width includes (twice) the pixel
      border. */
  float xdelt, ydelt;

  /** Rotation angle of the pixel grid around the optical axis going
      through the reference point [rad]. */
  float rota;

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

/** Return the x- and y-indices of the detector pixel corresponding to
    the specified x and y values. If the position is located outside
    the pixel area or on a pixel border, the returned indices are
    -1. The function also returns the relative position within the
    pixel. The latter values are required for the calculation of split
    events. */
void getGenDetAffectedPixel(const GenPixGrid* const grid, 
			    const double x, 
			    const double y,
			    int* const xi, 
			    int* const yi,
			    double* const xp, 
			    double* const yp);


#endif /* GENPIXGRID_H */
