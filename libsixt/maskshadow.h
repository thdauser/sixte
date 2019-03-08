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


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

#ifndef MASKSHADOW_H
#define MASKSHADOW_H 1

#include "sixt.h"
#include "vector.h"
#include "reconstruction.h"
#include "find_position.h"
#include "sourceimage.h"
#include "fft_array.h"
#include "eventarray.h"

////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////
typedef struct {
  double** map;    /**whole map of re-pixeled mask derived from ReconArray (values betw 0...1)*/
  double** shadow; /**part of re-pixeled mask corresponding to current source */
} MaskShadow;

/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////

MaskShadow* getEmptyMaskShadowElement(int* const status);

MaskShadow* getMaskShadowElement(int const Size1_map, int const Size2_map, int const Size1_shadow,
				 int const Size2_shadow, int* const status);

void getMaskRepix(ReconArray* recon, MaskShadow* ms, int shift);

void getMaskShadow(MaskShadow* ms, PixPositionList* ppl, SourceImage* sky_pixels, SquarePixels* det_pix,
		    ReconArray* r, const Vector* const nx, const Vector* const ny, double const distance);
void getMaskShadow2(MaskShadow* ms, struct wcsprm* wcs2, PixPositionList* ppl, float sky_crpix1, float sky_crpix2,
		    SquarePixels* det_pix, int const Size1, int const Size2, int const shift, int* const status);

double getNormalization1(MaskShadow* ms, ReadEvent* ea, SquarePixels* det_pix, int const xdiff, int const ydiff);
double getNormalization2(MaskShadow* ms, ReadEvent* ea, SquarePixels* det_pix, int const xdiff, int const ydiff);

void FreeMaskShadow(MaskShadow* ms,int const Size1);

#endif /* MASKSHADOW_H */
