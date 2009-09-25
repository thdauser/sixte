#ifndef ARCPIXELS_H
#define ARCPIXELS_H 1

#include "sixt.h"


/** One single ArcPixel. */
typedef struct {
  float charge; /**< Charge stored in this detector pixel. */
  double last_impact; /**< Time of the last photon impact in this pixel. */
} ArcPixel;


/** Array of ArcPixels. */
typedef struct {
  ArcPixel* array;
  
  /** Number of pixel rings. */
  int nrings; 
  /** Number of pixels in each individual ring. 
   * The array has nrings elements. */
  int* npixels;
} ArcPixels;


struct ArcPixelsParameters {
  int nrings;
  int* npixels;
};


/////////////////////////////////////////////////////////////////////


/** Initialization routine for the ArcPixels data structure. 
 * Sets the basic properties and allocates memory for the pixel array. */
int initArcPixels(ArcPixels*, struct ArcPixelsParameters*);


#endif /* ARCPIXELS_H */

