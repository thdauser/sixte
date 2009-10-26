#ifndef CODEDMASK_H
#define CODEDMASK_H 1

#include "sixt.h"
#include "photon.h"
#include "telescope.h"
#include "point.h"
#include "vector.h"


#define TRANSPARENT 1 /**< Pixel is transparent for X-rays. */
#define OPAQUE 0 /**< Pixel is opaque for X-rays. */


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

  /** Transparency of the mask.
   * The transparency is defined as the ratio #(transparent pixels)/#(all pixels). */
  double transparency;

  int naxis1, naxis2;    /**< Width of the image [pixel]. */
  double cdelt1, cdelt2; /**< Width of one pixel [m]. */
  double crpix1, crpix2; /**< Reference pixel [pixel] */
  double crval1, crval2; /**< Value at reference pixel [m] */

} CodedMask;


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////


/** Basic Constructor.
 * This constructor generates an empty CodedMask object.
 * It does not allocate any memory for the mask data. This has
 * to be done separately. */
CodedMask* getCodedMask(int* status);

/** Advanced Constructor.
 * The constructor obtains a new CodedMask object from the basic constructor
 * loads a coded mask from the specified file, and stores it in the newly
 * created object. */
CodedMask* getCodedMaskFromFile(char* filename, int* status);

/** Destructor. */
void freeCodedMask(CodedMask* mask);

/** Determine the impact position of a photon on the detector plane according to
 * the coded mask aperture.
 * The function return value is 1 if the photon passes the mask. If the photon is
 * absorbed by an opaque pixel in the mask, the return value is zero. */
int getCodedMaskImpactPos(struct Point2d* position, Photon* photon, CodedMask* mask, 
			  struct Telescope* telescope);


#endif /* CODEDMASK_H */
