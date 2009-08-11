#ifndef HEXAGONAL_PIXELS_H
#define HEXAGONAL_PIXELS_H 1

#include "sixt.h"
#include "point.h"
#include "genericdetector.h"


#define HTRS_N_PIXELS (37) // Total number of pixels in the HTRS array
#define HEXAGONAL_PIXELS_LINE_INDEX_OFFSET (7) // Offset of the line indices defining the pixels


/** One hexagonal detector pixel. */
typedef struct {
  float charge; /**< Charge stored in this detector pixel. */
} HexagonalPixel;


/** Different orientations for the 6 sub-triangles a hexagonal pixel consists of. */
typedef enum {
  TRIANGLE_ORIENTATION_RIGHT      =0,
  TRIANGLE_ORIENTATION_UPPER_RIGHT=1,
  TRIANGLE_ORIENTATION_UPPER_LEFT =2,
  TRIANGLE_ORIENTATION_LEFT       =3,
  TRIANGLE_ORIENTATION_LOWER_LEFT =4,
  TRIANGLE_ORIENTATION_LOWER_RIGHT=5
} SubTriangleDirection;


/** Contains information about the sub-triangle of a hexagonal pixel. */
typedef struct {
  /** Index of the hexagonal pixel. */
  int pixelindex;
  /** Orientation of the sub-triangle. */
  SubTriangleDirection orientation; 
} HexagonalPixelLineIndexInformation;


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

  /** Auxiliary array used to convert line indices to pixel indices. 
   * Addtionally to the pixel index corresponding to the given combination
   * of 3 line indices the array entry also provides information about the 
   * orientation of the corresponding sub-triangle of the hexagonal pixel. */
  HexagonalPixelLineIndexInformation lineIndexInformation
  [2*HEXAGONAL_PIXELS_LINE_INDEX_OFFSET +1]
  [2*HEXAGONAL_PIXELS_LINE_INDEX_OFFSET +1]
  [2*HEXAGONAL_PIXELS_LINE_INDEX_OFFSET +1];

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
int initHexagonalPixels(HexagonalPixels*, struct HexagonalPixelsParameters*);

/** Clean up the HexagonalPixels data structure. E.g. release allocated memory. */
void cleanupHexagonalPixels(HexagonalPixels* sp);

/** Clear the array of HexagonalPixels. */
inline void clearHexagonalPixels(HexagonalPixels*);

/** Determine the Pixel index of the pixel(s) that is (are) affected by the photon
 * impact at the given position. 
 * The routine determines the primary impact pixel and, if appropriate, also the 
 * split partner in a double split event, if the selected charge cloud size is none-zero,
 * and the event is close enough to the pixel edge. 
 * The accroding charge fractions are also determined by the routine.
 * The generated charge is distributed among the two pixels according to a Gaussian shape
 * charge cloud model. 
 * The return valud of the function is the number of split partners, i.e., it is either
 * 1 or 2. */
int getHexagonalPixelSplits(HexagonalPixels* hp, GenericDetector* gd, 
			    struct Point2d position, int* pixel, double* fraction);


#endif /* HEXAGONAL_PIXELS_H */

