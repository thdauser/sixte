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

#include "fft_array.h"

fftw_complex* FFTOfArray_1d(double* Image1d, int ImageSize1, int ImageSize2, int type)
{
  fftw_complex* Input;
  fftw_complex* Output;
  fftw_plan plan;
  int ii,jj;

  //Memory-Allocation
  Input=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ImageSize1*ImageSize2);
  Output=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ImageSize1*ImageSize2);

  //Copy Image1d of type double to Input-Array of type fftw_complex

  for(ii=0; ii<ImageSize1; ii++){
    for(jj=0; jj<ImageSize2; jj++){
    Input[ii+ImageSize1*jj][0]=Image1d[ii+ImageSize1*jj];
    Input[ii+ImageSize1*jj][1]=0.;
    }
  }

  plan=fftw_plan_dft_2d(ImageSize1, ImageSize2, Input, Output, type, FFTW_ESTIMATE);

  fftw_execute(plan);

  fftw_destroy_plan(plan);
  fftw_free(Input);

  return(Output);
}


fftw_complex* FFTOfArray(fftw_complex* Input, int ImageSize1, int ImageSize2, int type)
{
  fftw_complex* Output;
  fftw_plan plan;

  //Memory-Allocation
  Output=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ImageSize1*ImageSize2);

  plan=fftw_plan_dft_2d(ImageSize1, ImageSize2, Input, Output, type, FFTW_ESTIMATE);

  fftw_execute(plan);

  fftw_destroy_plan(plan);
  fftw_free(Input);

  return(Output);
}
