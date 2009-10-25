#ifndef CODEDMASK_H
#define CODEDMASK_H 1

#include "sixt.h"
#include "photon.h"


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** Contains a coded mask. */
typedef struct {
  /** The actual coded mask data.
   * This array respresents the individual pixels of the Coded Mask.
   * For transparent pixels the array value is 1, for opaque pixels it
   * is zero. */
  int** map;

  /** A 1-dimensional list of the transparent pixels. 
   * This list contains the transparent pixels only.
   * For each transparent pixel the indices x and y in the 2-dimensional
   * map are stored in this list. 
   * The list is used to determine randomly a pixel where the photon is
   * transmitted. */
  int** transparent_pixels;
  /** Number of transparent pixels. */
  int n_transparent_pixels;

  /** Opacity of the mask.
   * The opacity is defined as the ratio #(opaque pixels)/#(all pixels). */
  double opacity;

  int naxis1, naxis2;    /**< Width of the image [pixel]. */
  double cdelt1, cdelt2; /**< Width of one pixel [m]. */
  double crpix1, crpix2; /**< [pixel] */
  double crval1, crval2; /**< [m] */
} CodedMask;


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////


/** Basic Constructor.
 * This constructor generates an empty CodedMask object.
 * It does not allocate any memory for the mask data. This has
 * to be done separately. */
CodedMask* getCodedMask(int* status);

/** Destructor. */
void freeCodedMask(CodedMask* mask);


#endif /* CODEDMASK_H */
