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

#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H 1


#include "codedmask.h"
#include "sixt.h"
#include "squarepixels.h"
#include "repix.h"

/////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////
typedef struct {
  double** Rmap; //the reconstruction-array data built from mask-array-data.
  double open_fraction; //Transparency of the mask.
                        //Ratio #(transparent pixels)/#(all pixels)
  int naxis1, naxis2;   // Width of the image [pixel]
}ReconArray;


/////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////
ReconArray* newReconArray(int* const status);

ReconArray* getReconArray(CodedMask* mask,int type, SquarePixels* detector_pixels, int* const status);

double* SaveReconArray1d(ReconArray* recon, int* status);

void FreeReconArray1d(double* ReconArray1d);

void FreeReconArray(ReconArray** const recon);

//void FreeReconImage(ReconArray* recon);


#endif
