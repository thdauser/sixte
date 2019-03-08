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
   Copyright 2018-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "tescrosstalk.h"


/** Assign resonance frequencies to advdet pixels */
void get_resonance_frequencies(AdvDet* det, int* status){
  // this code uses a lot of the logic from get_readout_channels in crosstalk.c

  // first, check whether a file was given (resfreq is initialized to freq)
  // if not, it would be an empty string
  // an empty string starts with the terminator; this is faster than strcmp
  if (det->channel_resfreq_file[0] == '\0'){
    // don't need to do anything
    return;
  }
  // get the complete filename
  char fullfilename[MAXFILENAME];
  strcpy(fullfilename, det->filepath);
  strcat(fullfilename, det->channel_resfreq_file);

  // get the channel list
  channel_list* chans = load_channel_list(fullfilename, status);
  CHECK_STATUS_VOID(*status);

  // check if the file agrees with the number of pixels
  if (chans->len != det->npix){
    printf("*** error: number of pixels from channel_file %s (%i) is not equal total number of pixels (%i)!\n",
           fullfilename,chans->len,det->npix);
    *status=EXIT_FAILURE;
    return;
  }

  for (int ii=0; ii < chans->len; ii++){
    // check if PixID makes sense (PixID starts at 0)
    if ( (chans->pixid[ii]<0) || (chans->pixid[ii] > det->npix-1)){
      printf("*** error: Pixel-ID %i does not belong to the detector \n", chans->pixid[ii]);
      *status = EXIT_FAILURE;
      return;
    }

    // get the relevant pixel
    AdvPix* pix = &(det->pix[chans->pixid[ii]]);

    // set resfreq of the pixel
    pix->resfreq = chans->freq[ii];
    if (pix->freq <= 0.0){
      printf("*** warning: assigning unrealistic frequency (%.3e) to pixel %i\n", pix->freq, pix->pindex);
    }
  }
  // free channel list
  free_channel_list(&chans);
}


/** Alloc an empty FDM system */
FDMSystem* newFDMSystem(int num_pixels, int* status){

  FDMSystem* fdmsys = (FDMSystem*) malloc(sizeof (FDMSystem));
  CHECK_MALLOC_RET_NULL_STATUS(fdmsys, *status);

  fdmsys->num_pixels = num_pixels;

  fdmsys->omega_array = (double*) malloc(num_pixels*sizeof(double));
  CHECK_MALLOC_RET_NULL_STATUS(fdmsys->omega_array, *status);

  fdmsys->res_omega_array = (double*) malloc(num_pixels*sizeof(double));
  CHECK_MALLOC_RET_NULL_STATUS(fdmsys->res_omega_array, *status);

  // set up the Z_array
  fdmsys->Z_array = (double **)malloc(num_pixels* sizeof(double));
  CHECK_MALLOC_RET_NULL_STATUS(fdmsys->Z_array, *status);
  int ii;
  for (ii=0; ii<num_pixels; ii++){
    fdmsys->Z_array[ii] = (double *)malloc(num_pixels*sizeof(double));
    CHECK_MALLOC_RET_NULL_STATUS(fdmsys->Z_array[ii], *status);
  }

  fdmsys->L_Common = 0.;
  fdmsys->C_Common = -1.;

  fdmsys->X_L = (double*) malloc(num_pixels*sizeof(double));
  CHECK_MALLOC_RET_NULL_STATUS(fdmsys->X_L, *status);

  fdmsys->capFac = (double*) malloc(num_pixels*sizeof(double));
  CHECK_MALLOC_RET_NULL_STATUS(fdmsys->capFac, *status);

  return fdmsys;

}


/** Initialize an FDM system */
void init_FDMSystem(Channel* chan, double L_Common, double C_Common, double TTR, int* status){

  // get an empty FDMSystem
  FDMSystem* sys = newFDMSystem(chan->num_pixels, status);

  int npix = chan->num_pixels; // we'll need this a lot
  sys->L_Common = L_Common;
  sys->C_Common = C_Common;
  int ii, jj; // for looping

  // assign the frequency array and X_L
  for (ii=0;ii<npix;ii++){
    sys->omega_array[ii] = 2*M_PI * chan->pixels[ii]->freq;
    sys->res_omega_array[ii] = 2*M_PI * chan->pixels[ii]->resfreq;
    sys->X_L[ii] = sys->omega_array[ii] * sys->L_Common / TTR /TTR; // down-transformed by TTR
  }

  // complex impedances in Z_array
  // index structure: Z_i^{\omega_j} = Z_array[j][i]
  for (ii=0;ii<npix;ii++){
    for (jj=0;jj<npix;jj++){
      // need to use resonance frequencies for pixel ii
      sys->Z_array[jj][ii] = (chan->pixels[ii]->tes->Lfilter + L_Common) * (sys->omega_array[jj] * sys->omega_array[jj] - sys->res_omega_array[ii]*sys->res_omega_array[ii]) / sys->omega_array[jj] /TTR /TTR; // down-transformed by TTR
    }
  }

  // include common capacitance
  // capfac[i_freq]!
  if (npix != 2 || sys->C_Common <= 0.) {
    // all capFacs are zero
    for (ii=0;ii<npix;ii++){
      sys->capFac[ii] = 0.;
    }
  }
  if (sys->C_Common > 0.){
    if (npix != 2){
      headas_chat(0,"Common Capacitance Crosstalk is only supported for two pixels!\n");
    } else {
      for (ii=0;ii<npix;ii++){
        for (jj=0; jj<npix;jj++){
          if (ii==jj) continue;
          double X_C_inv = -1.*(sys->C_Common*sys->omega_array[jj])*TTR*TTR; //X_C^-1. UP-transformed, since it's inverse
          sys->capFac[jj] = (sys->Z_array[jj][ii] - 2 * sys->X_L[jj])*(X_C_inv);
        }
      }
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
        gsl_complex denominator = gsl_complex_rect((chan->pixels[i_pix]->tes->RT+chan->pixels[i_pix]->tes->Reff)*(1.+chan->fdmsys->capFac[i_freq]), chan->fdmsys->Z_array[i_freq][i_pix] - chan->fdmsys->X_L[i_freq]);
        gsl_complex numerator = gsl_complex_rect((chan->pixels[i_freq]->tes->RT+chan->pixels[i_freq]->tes->Reff)*(1.+chan->fdmsys->capFac[i_freq]), chan->fdmsys->Z_array[i_freq][i_freq] - chan->fdmsys->X_L[i_freq]);
        Icalc = gsl_complex_mul_real(gsl_complex_div(numerator,denominator), chan->pixels[i_freq]->tes->I0);
      }
      // write the calculated current into the corresponding tes pixels
      // carrier overlap: add the current to Ioverlap of that frequency
      chan->pixels[i_freq]->tes->Ioverlap = gsl_complex_add(chan->pixels[i_freq]->tes->Ioverlap,Icalc);
      // common impedance: add abs(I)^2 to Pcommon of the current pixel
      chan->pixels[i_pix]->tes->Pcommon += gsl_complex_abs2(Icalc);
      //chan->pixels[i_pix]->tes->Pcommon += GSL_REAL(Icalc)*chan->pixels[i_pix]->tes->V0;
    }
  }
  // finish Common Impedance power and find phase rotation
  for (i_pix=0; i_pix<chan->num_pixels; i_pix++){
    // P = I^2 *R
    chan->pixels[i_pix]->tes->Pcommon *= chan->pixels[i_pix]->tes->RT;

    // calculate bias voltage from currents
    // contribution from the on-resonance pixel
    gsl_complex Vbias = gsl_complex_mul_real(gsl_complex_rect(chan->pixels[i_pix]->tes->RT+chan->pixels[i_pix]->tes->Reff , chan->fdmsys->Z_array[i_pix][i_pix]) , chan->pixels[i_pix]->tes->I0);
    // contribution from the other pixels
    Vbias = gsl_complex_add(Vbias, gsl_complex_mul_imag(chan->pixels[i_pix]->tes->Ioverlap,chan->fdmsys->X_L[i_pix]));
    // to get the actual I and q channel, all currents need to be rotated by the angle of V
    chan->pixels[i_pix]->tes->theta_Vb = gsl_complex_arg(Vbias);
  }
}
