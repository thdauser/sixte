/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, Mirjam Oertel, FAU
*/

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
  /** The actual coded mask data. This array respresents the
      individual pixels of the Coded Mask. For transparent pixels the
      array value is 1, for opaque pixels it is zero. */
  int** map;

  /** A 1-dimensional list of the transparent pixels. This list
      contains the transparent pixels only. For each transparent pixel
      the indices x and y in the 2-dimensional map are stored in this
      list. The list is used to determine randomly a pixel where the
      photon is transmitted. */
  int** transparent_pixels;
  /** Number of transparent pixels. */
  int n_transparent_pixels;

  /** Transparency of the mask. The transparency is defined as the
      ratio #(transparent pixels)/#(all pixels). */
  double transparency;

  int naxis1, naxis2;    /**< Width of the image [pixel]. */
  double cdelt1, cdelt2; /**< Width of one pixel [m]. */
  double crpix1, crpix2; /**< Reference pixel [pixel] */
  double crval1, crval2; /**< Value at reference pixel [m] */

} CodedMask;


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////


/** Basic Constructor. This constructor generates an empty CodedMask
    object. It does not allocate any memory for the mask data. This
    has to be done separately. */
CodedMask* getCodedMask(int* const status);

/** Advanced Constructor. The constructor obtains a new CodedMask
    object from the basic constructor loads a coded mask from the
    specified file, and stores it in the newly q created object. */
CodedMask* getCodedMaskFromFile(const char* const filename, int* const status);

/** Destructor. */
void destroyCodedMask(CodedMask** const mask);


/** Determine the impact position of a photon in the detector plane.
    The function return value is 1 if the photon passes an transparent mask-pixel
    and doesn't hit the walls. Otherwise the return value is zero. */

int getImpactPos (struct Point2d* const position,
		  const Vector* const phodir,
		  const CodedMask* const mask, 
		  const struct Telescope* const telescope,
		  const Vector* const nz,
		  const float distance,
		  const float x_det,
		  const float y_det,
		  int* const status);

int getImpactPos_wcs (struct wcsprm* wcs, struct Point2d* const position, const CodedMask* const mask,
		   double const photon_ra, double const photon_dec, float const det_pixelwidth,
		   const float x_det, const float y_det, int* const status);

int getImpactPos_protoMirax(struct wcsprm* wcs, struct wcsprm* wcs2, struct Point2d* const position,
			   const CodedMask* const mask, double const photon_ra, double const photon_dec,
			   float const det_pixelwidth, const float det_width, const float x_det,
			    const float y_det, const float wall, int* const status);


#endif /* CODEDMASK_H */
