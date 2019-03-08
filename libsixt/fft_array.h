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

#ifndef FFT_ARRAY_H
#define FFT_ARRAY_H1

#include "sixt.h"
#include "fftw3.h"

////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////

//performs a fft of an input array which has to be in row-major format (1d-image)
//type +1 equals FFTW_BAKWARD;type -1 equals FFTW_FORWARD
fftw_complex* FFTOfArray_1d(double* Image1d, int ImageSize1, int ImageSize2, int type);

//performs a fft of an input array which is of type fftw_complex
//type +1 equals FFTW_BAKWARD;type -1 equals FFTW_FORWARD
fftw_complex* FFTOfArray(fftw_complex* Input, int ImageSize1, int ImageSize2, int type);

#endif
