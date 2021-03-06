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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef TESTIMG_H
#define TESTIMG_H 1

#include "sixt.h"
#include "eventarray.h"
#include "maskshadow.h"
#include "fft_array.h"
#include "projectedmask.h"
#include "skyimage.h"


#define TMASKSHADOW 1     //MaskShadow* ms->shadow
#define TREADEVENT 2      //ReadEvent* ea->EventArray
#define TMASKSHADOWMAP 3  //MaskShadow* ms->map
#define TPROJMASK 4       //ProjectedMask* pm->map
#define TMASKMAP 5        //CodedMask* cm->map
#define TSOURCEIMG 6      //SourceImage* si->pixel
#define TSKYIMG 7         //SkyImage* si->pixel OR si->ra OR si->dec

////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
// Function Declarations.
////////////////////////////////////////////////////////////////////////

//if type==TREADEVENT -> shiftX=ea->naxis1/2+xdiff -> sizeX=det_pix->xwidth
//else -> shiftX=0
void createTestImg(void* arg, int type, int sizeX, int sizeY, int shiftX,
		   int shiftY, char* filename, int* const status);


void testFitsImage1d(double* Image1d, char* filename, int Size1, int Size2,
		     int* const status);

//type=0 for real-part, type=1 for imaginary part
void createTestImg_fft_part(fftw_complex* Image, char* filename, int type, int Size1, int Size2, int* const status);

#endif /* TESTIMG_H */
