#ifndef SQUARE_PIXELS_H
#define SQUARE_PIXELS_H 1

#include "sixt.h"
#include "point.h"
#include "genericdetector.h"


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** One square detector pixel. */
typedef struct {
  float charge; /**< Charge stored in this detector pixel. */
} SquarePixel;


/** Data specific for detectos with an array of square pixels. */
typedef struct {
  int xwidth; /**< Width of the detector pixel array in x-direction 
	       * (number of [integer pixels]). */
  int ywidth; /**< Width of the detector pixel array in y-direction
	       * (number of [integer pixels]). */

  int xoffset; /**< Offset of the detector pixel array in x-direction in [integer pixels].
	       * The physical origin of the detector (at the center of the detector) 
	       * has the array x-index 'xoffset'. */
  int yoffset; /**< Offset of the detector pixel array in y-direction in [integer pixels].
	       * The physical origin of the detector (at the center of the detector) 
	       * has the array y-index 'yoffset'. */

  double xpixelwidth; /**< Width of a single pixel in x-direction [m]. */
  double ypixelwidth; /**< Width of a single pixel in y-direction [m]. */

  SquarePixel** array; /**< Array of square pixels. */

} SquarePixels;


struct SquarePixelsParameters {
  int xwidth, ywidth;
  double xpixelwidth, ypixelwidth;
};


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////


/** Initialization routine for the SquarePixels data structure. 
 * Sets the basic properties and allocates memory for the pixel array. */
int initSquarePixels(SquarePixels*, struct SquarePixelsParameters*);

/** Clean up the SquarePixels data structure. E.g. release allocated memory. */
void cleanupSquarePixels(SquarePixels* sp);

/** Clear the array of SquarePixels. */
inline void clearSquarePixels(SquarePixels*);

/** Clear on line of Pixels. A line is defined to have constant detector x-coordinate. */
inline void clearLineSquarePixels(SquarePixels* sp, const int line);

/** Determine the split ratios of a photon impact on an array of square pixels. */
int getSquarePixelsSplits(SquarePixels*, GenericDetector*, struct Point2d position, 
			  int* x, int* y, double* fraction);


#endif /* SQUARE_PIXELS_H */

