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
int initArcPixels(ArcPixels*, struct ArcPixelsParameters*);

/** Clean up the ArcPixels data structure. E.g. release allocated memory. */
void cleanupArcPixels(ArcPixels*);

/** Clear the array of ArcPixels. */
inline void clearArcPixels(ArcPixels*);

/** Determine the ArcPixel that contains the specified position. 
 * Return values are the detector ring and the pixel index within this ring. */
void getArcPixel(ArcPixels*, struct Point2d position, int* ring, int* pixel);


#endif /* ARCPIXELS_H */

