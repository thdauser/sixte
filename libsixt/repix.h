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

#ifndef REPIX_H
#define REPIX_H 1

#include "sixt.h"
#include "eventarray.h"
#include "reconstruction.h"
#include "projectedmask.h"

#define TMASKSHADOW 1     //MaskShadow* ms->shadow
#define TREADEVENT 2      //ReadEvent* ea->EventArray
#define TMASKSHADOWMAP 3  //MaskShadow* ms->map
#define TPROJMASK 4       //ProjectedMask* pm->map
#define TMASKMAP 5        //CodedMask* cm->map,ReconArray* ra->Rmap


/////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////

//re-pixels from bigger pixels to smaller ones that fit without reminder into former bigger ones
//only for square pixels!

void repixNoReminder(void* arg_dataBeforeRepix, void* arg_dataAfterRepix, int type,
		     int Size1, int Size2, double pixelwidth_big, double pixelwidth_small);

  //arg_dataBeforeRepix: pointer to array which shall be re-pixeled
  //arg_dataAfterRepix: pointer to array where re-pixeled version is stored
  //type: data type of both arrays
  //Size1/2: sizes of array which shall be re-pixeled (in pixels)
  //pixelwidth_big: size of pixels of given array before re-pixelization (in meters)
  //pixelwidth_small: size of pixels to which the array shall be re-pixeled to (in meters)


//re-pixels from bigger pixels to smaller ones that fit not without reminder into former bigger ones
//the re-pixeled array is hence blurred due to approximation of the data

void repixWithReminder(void* arg_dataBeforeRepix, void* arg_dataAfterRepix, int type,
		       int Size1, int Size2, double pixelwidth_BR_big, double pixelwidth_BR_small,
		       double pixelwidth_small, double MinVal);
  //pixelwidth_BR_big/small: width of former (BeforeRepix) big pixels; in case of varying pixelsize in
    //input array(even: big, odd: small(but still bigger than RePixSize))give both values, else
    //give the same value for both parameter
  //MinVal: min value that may occur in the re-pixeled array; the whole raw new array (arg_dataAfterRepix)
    //must be initialized to that value before handling it to the repix-fct

double getRepixValue(SquarePixels* det_pix, double pixelsize1, /*double pixelsize2,*/ int factor);


#endif /* REPIX_H */
