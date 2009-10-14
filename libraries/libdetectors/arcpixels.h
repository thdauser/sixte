#ifndef ARCPIXELS_H
#define ARCPIXELS_H 1

#include "sixt.h"
#include "point.h"
#include "genericdetector.h"


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** One single ArcPixel. */
typedef struct {
  float charge; /**< Charge stored in this detector pixel. */
  double last_impact; /**< Time of the last photon impact in this pixel. */
} ArcPixel;


/** Array of ArcPixels. */
typedef struct {
  ArcPixel** array;
  
  /** Number of pixel rings. */
  int nrings; 
  /** Number of pixels in each individual ring. 
   * The array has nrings elements. */
  int* npixels;

  /** Radii of the individual detector rings. 
   * The array has nrings elements each specifying the outer radius of the
   * respective pixel ring. The inner radius of the ring i is given by the radius of the 
   * ring i-1 (for i>0). */
  double* radius;
  /** Offset angle of the first pixel in each ring ([rad]).
   * In general the pixels may not start the 3 o'clock position in the rings, but
   * the border of pixel with the index 0 may have an angular offset from the 3 o'clock
   * position. These offsets are measured in counter-clock-wise direction. */
  double* offset_angle;

} ArcPixels;


struct ArcPixelsParameters {
  int nrings;
  int* npixels;
  double* radius;
  double* offset_angle;
};


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////


/** Initialization routine for the ArcPixels data structure. 
 * Sets the basic properties and allocates memory for the pixel array. */
int initArcPixels(ArcPixels* ap, struct ArcPixelsParameters* app);

/** Clean up the ArcPixels data structure. E.g. release allocated memory. */
void cleanupArcPixels(ArcPixels* ap);

/** Clear the array of ArcPixels. */
inline void clearArcPixels(ArcPixels* ap);

/** Determine the ArcPixel that contains the specified position. 
 * The return value is the absolute pixel index (i.e., within the global
 * numbering of all pixels of this detector) of the affected pixel. */
void getArcPixel(ArcPixels* ap, struct Point2d position, int* pixel);

/** Determine the absolute pixel index from a given ring and the pixel 
 * number within this ring. */
int getArcPixelIndex(ArcPixels* ap, int ring, int number);

#endif /* ARCPIXELS_H */

