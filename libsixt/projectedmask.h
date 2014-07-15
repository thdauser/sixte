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

#ifndef PROJECTEDMASK_H
#define PROJECTEDMASK_H 1

#include "sixt.h"
#include "codedmask.h"
#include "squarepixels.h"

////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////
typedef struct{
  double** map;  /* Map of the projected mask with different alternating pixel-sizes.
		    The pixels of size 'pixelwidth1 x pixelwidth1' contain the former
		    mask values (1 equals tranparent, 0 equals opaque).
		    The pixels of size 'pixelwidth2 x pixelwidth1' contain the intermediate values
		    given by the two surrounding 'real' former mask pixels and the pixels of size 
		    'pixelwidth2 x pixelwidth2' contain the intermediate values of the four
		    surrounding 'intermediate pixels'. */

  double pixelwidth1;  /* Projected width with shift in mask and det plane.
			  proj_dist*(mask_pixelwidth-det_pixelwidth).
		          In each line: all odd pixels.*/
  double pixelwidth2;  /* Projected width with shift only in mask plane.
			  proj_dist*det_pixelwidth.
		          In each line: all even pixels, each line starts and ends
		          with such a pixel.*/

  int naxis1;  /*  Width of the image [pixel]. 2*amount_of_mask_pixels+1 */
  int naxis2;  /*  Width of the image [pixel]. 2*amount_of_mask_pixels+1 */

} ProjectedMask;

/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////

ProjectedMask* newProjectedMask(int* const status);

ProjectedMask* getProjectedMask(const CodedMask* const mask, SquarePixels* det_pix, 
				const double proj_dist, int* const status);

#endif /* PROJECTEDMASK_H */
