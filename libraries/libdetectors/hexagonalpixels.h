#ifndef HEXAGONAL_PIXELS_H
#define HEXAGONAL_PIXELS_H 1

#include "sixt.h"
#include "point.h"
#include "genericdetector.h"


/** One hexagonal detector pixel. */
typedef struct {
  float charge; /**< Charge stored in this detector pixel. */
} HexagonalPixel;


/** Data specific for detectos with an array of hexagonal pixels. */
typedef struct {
  int npixels; /**< Total number of pixels. Currently has to be 37. */

  /** Length of one egde of a hexagonal pixel (in [m]). 
   * The length of a pixel edge 'a' is related to the height of the pixels 'h'. */
  double a; 

  /** Height of one of the six equilateral triangles that define
   * the hexagonal structure of the pixels (in [m]). 
   * The value of 'h' is related to the length of the pixel egdes 'a'. */
  double h; 

  HexagonalPixel* array; /**< 1-dimensional array of hexagonal pixels. */

} HexagonalPixels;


// Dimensions of a hexagonal pixel:
//
//       -------        -
//      /   |   \       |
//     /    |    \      |
//    /     |     \     |
//    \     |2h   /     | pixelwidth
//     \    |    /      |
//      \   |   /       |
//       -------        -
//          a
// 
// with:   a = 2 * h / sqrt(3)
// and :   h = pixelwidth / 2


struct HexagonalPixelsParameters {
  int npixels;
  double pixelwidth;
};


/////////////////////////////////////////////////////////////////////


/** Initialization routine for the HexagonalPixels data structure. 
 * Sets the basic properties and allocates memory for the pixel array. */
//int initHexagonalPixels(HexagonalPixels*, struct HexagonalPixelsParameters*);

/** Clean up the HexagonalPixels data structure. E.g. release allocated memory. */
//void cleanupHexagonalPixels(HexagonalPixels* sp);

/** Clear the array of HexagonalPixels. */
//inline void clearHexagonalPixels(HexagonalPixels*);

/** Clear on line of Pixels. A line is defined to have constant detector x-coordinate. */
//inline void clearLineHexagonalPixels(HexagonalPixels* sp, const int line);

/** Determine the split ratios of a photon impact on an array of hexagonal pixels. */
//int getHexagonalPixelsSplits(HexagonalPixels*, GenericDetector*, struct Point2d position, 
//			     int* x, int* y, double* fraction);


#endif /* HEXAGONAL_PIXELS_H */

