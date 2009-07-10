#ifndef SQUARE_PIXELS_H
#define SQUARE_PIXELS_H 1

#include "sixt.h"

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


/** Initialization routine for the SquarePixels data structure. 
 * Sets the basic properties and allocates memory for the pixel array. */
int initSquarePixels(SquarePixels*, struct SquarePixelsParameters*);


#endif /* SQUARE_PIXELS_H */

