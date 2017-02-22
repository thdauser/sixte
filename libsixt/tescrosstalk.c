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


   Copyright 2017 Christian Kirsch, ECAP
*/

#include "tescrosstalk.h"


/** Alloc an empty FDM system */
FDMSystem* newFDMSystem(int num_pixels, int* status){
  
  FDMSystem* fdmsys = (FDMSystem*) malloc(sizeof (FDMSystem));
  CHECK_MALLOC_RET_NULL_STATUS(fdmsys, *status);

  fdmsys->num_pixels = num_pixels;

  fdmsys->omega_array = (double*) malloc(num_pixels*sizeof(double));
  CHECK_MALLOC_RET_NULL_STATUS(fdmsys->omega_array, *status);

  // set up the Z_array 
  fdmsys->Z_array = (double **)malloc(num_pixels* sizeof(double));
  CHECK_MALLOC_RET_NULL_STATUS(fdmsys->Z_array, *status);
  int ii;
  for (ii=0; ii<num_pixels; ii++){
    fdmsys->Z_array[ii] = (double *)malloc(num_pixels*sizeof(double));
    CHECK_MALLOC_RET_NULL_STATUS(fdmsys->Z_array[ii], *status);
  }

  fdmsys->L_Common = 0;

  fdmsys->u_LC = (double*) malloc(num_pixels*sizeof(double));
  CHECK_MALLOC_RET_NULL_STATUS(fdmsys->u_LC, *status);

  return fdmsys;

}


/** Initialize an FDM system */
void init_FDMSystem(Channel* chan, double L_Common, double TTR, int* status){
  
  // get an empty FDMSystem
  FDMSystem* sys = newFDMSystem(chan->num_pixels, status);

  int npix = chan->num_pixels; // we'll need this a lot
  sys->L_Common = L_Common;
  int ii, jj; // for looping

  // assign the frequency array and u_LC
  for (ii=0;ii<npix;ii++){
    sys->omega_array[ii] = 2*M_PI * chan->pixels[ii]->freq;
    sys->u_LC[ii] = sys->omega_array[ii] * sys->L_Common / TTR /TTR; // down-transformed by TTR
  }

  // complex impedances in Z_array
  // index structure: Z_i^{\omega_j} = Z_array[j][i]
  for (ii=0;ii<npix;ii++){
    for (jj=0;jj<npix;jj++){
      sys->Z_array[jj][ii] = (chan->pixels[ii]->tes->Lfilter + L_Common) * (sys->omega_array[jj] * sys->omega_array[jj] - sys->omega_array[ii]*sys->omega_array[ii]) / sys->omega_array[jj] /TTR /TTR; // down-transformed by TTR
    }
  }
  // assign the FDMSystem to the channel
  chan->fdmsys = sys;
}

/** Solve this channel's FDM system for the current channel state and output to the pixels */
void solve_FDM(Channel *chan){

  int i_freq, i_pix;
  // set Ioverlap and Pcommon to 0
  for (int ii=0; ii<chan->num_pixels; ii++){
    chan->pixels[ii]->tes->Ioverlap = gsl_complex_rect(0.,0.);
    chan->pixels[ii]->tes->Pcommon = 0;
  }
  // calculate off resonance current I_{i_pix}^{\omega_{i_freq}}
  gsl_complex Icalc;
  for (i_freq=0; i_freq<chan->num_pixels; i_freq++){
    for (i_pix=0; i_pix<chan->num_pixels; i_pix++){
      if (i_pix == i_freq) {
        continue; // only off-resonance currents
        //Icalc = gsl_complex_rect(0,0);
      } else { 
        gsl_complex denominator = gsl_complex_rect(chan->pixels[i_pix]->tes->RT, chan->fdmsys->Z_array[i_freq][i_pix] - chan->fdmsys->u_LC[i_freq]);
        gsl_complex numerator = gsl_complex_rect(chan->pixels[i_freq]->tes->RT, -1.*chan->fdmsys->u_LC[i_freq]);
        Icalc = gsl_complex_mul_real(gsl_complex_div(numerator,denominator), chan->pixels[i_freq]->tes->I0);
      }
      // write the calculated current into the corresponding tes pixels
      // carrier overlap: add the current to Ioverlap of that frequency
      chan->pixels[i_freq]->tes->Ioverlap = gsl_complex_add(chan->pixels[i_freq]->tes->Ioverlap,Icalc);
      // common impedance: add abs(I)^2 to Pcommon of the current pixel
      chan->pixels[i_pix]->tes->Pcommon += gsl_complex_abs2(Icalc);
    }
  }
  // P = I^2 *R
  for (i_pix=0; i_pix<chan->num_pixels; i_pix++){
    chan->pixels[i_pix]->tes->Pcommon *= chan->pixels[i_pix]->tes->RT;
  }
}
