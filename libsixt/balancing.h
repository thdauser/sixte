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

#ifndef BALANCING_H
#define BALANCING_H 1


#include "sixt.h"
#include "squarepixels.h"
#include "reconstruction.h"
#include "eventlist.h"


/////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////

typedef struct {
  double** Bmap; //the balancing-array data built from reconstruction-array-data.
  int naxis1, naxis2;    // Width of the image [pixel]
}BalancingArray;



/////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////

BalancingArray* newBalancingArray(int* const status);

BalancingArray* getBalancingArray(ReconArray* recon, SquarePixels* detector_pixels, EventList* ef, int* const status);

#endif
