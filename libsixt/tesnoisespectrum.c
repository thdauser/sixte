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


   Copyright 2014 Jelle de Plaa, SRON, Thorsten Brand, FAU 
*/


#include "tesnoisespectrum.h"


NoiseSpectrum* newNoiseSpectrum(AdvDet *det, 
				int* const status) {
    
    NoiseSpectrum* Noise=malloc(sizeof(NoiseSpectrum));
    if(Noise==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for noise spectrum failed");
      CHECK_STATUS_RET(*status, Noise);
    }
    int ii;
    
    /* Set noise parameters */
    Noise->WhiteRMS=det->TESNoise->WhiteRMS;  /* White Noise normalisation */
    Noise->H0=det->TESNoise->H0;              /* Noise filter normalisation */
    Noise->Nz=det->TESNoise->Nz;              /* Number of zeros in the noise spectrum */
    Noise->Np=det->TESNoise->Np;              /* Number of poles in the noise spectrum */

    Noise->Poles=NULL;
    Noise->Zeros=NULL;
    
    if(Noise->Np>0){
      Noise->Poles=(double*)malloc(Noise->Np*sizeof(double));
      if(Noise->Poles==NULL){
	*status=EXIT_FAILURE;
	SIXT_ERROR("memory allocation for noise spectrum failed");
	CHECK_STATUS_RET(*status, Noise);
      }
      for(ii=0; ii<Noise->Np; ii++){
	Noise->Poles[ii]=det->TESNoise->Poles[ii];
      }
    }
    if(Noise->Nz>0){
      Noise->Zeros=(double*)malloc(Noise->Nz*sizeof(double));
      if(Noise->Zeros==NULL){
	*status=EXIT_FAILURE;
	SIXT_ERROR("memory allocation for noise spectrum failed");
	CHECK_STATUS_RET(*status, Noise);
      }
      for(ii=0; ii<Noise->Nz; ii++){
	Noise->Zeros[ii]=det->TESNoise->Zeros[ii];
      }
    }
    return Noise;
}    


NoiseBuffer* newNoiseBuffer(int* const status, 
			    int *NumberOfPixels) 
{
    int i;
    
    /* Set Buffer properties */
    
    NoiseBuffer* NBuffer=(NoiseBuffer*)malloc(sizeof(NoiseBuffer));
    if(NBuffer==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for noise buffer failed");
      CHECK_STATUS_RET(*status, NBuffer);
    }
    
    NBuffer->BufferSize=NOISEBUFFERSIZE;
    NBuffer->NPixel=*NumberOfPixels;
    NBuffer->Buffer=(double**)malloc(NBuffer->BufferSize*sizeof(double*));
    if(NBuffer->Buffer==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for NBuffer Buffer failed");
      CHECK_STATUS_RET(*status, NBuffer);
    }
    for (i=0;i<NBuffer->BufferSize;i++) {
      NBuffer->Buffer[i]=(double*)malloc(NBuffer->NPixel*sizeof(double));
      if(NBuffer->Buffer[i]==NULL){
	*status=EXIT_FAILURE;
	SIXT_ERROR("memory allocation for NBuffer Buffer failed");
	CHECK_STATUS_RET(*status, NBuffer);
      }
    }
    
    return NBuffer;
}


int genNoiseSpectrum(NoiseSpectrum* Noise, 
		     NoiseBuffer* NBuffer, 
		     double *SampFreq, 
                     int* const status) 
{
    
    const gsl_rng_type * T;
    gsl_rng * r;
    double U, G, sigma;
    const double pi=3.14159265359;
    int i,j, k;
    
    fftw_complex *in, *out;
    fftw_plan p;
    fftw_complex Ze, Po, H; /* Products of Zeros & Poles */
    double w, f;               /* Omega and frequency*/
    
    /* gsl_rng_set(); */
    
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    sigma=1.;
    
    in=(fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NBuffer->BufferSize);
    if(in==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for fftw_complex failed");
      CHECK_STATUS_RET(*status, *status);
    }
    out=(fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NBuffer->BufferSize);
    if(out==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for fftw_complex failed");
      CHECK_STATUS_RET(*status, *status);
    }
    
    for (j=0; j<NBuffer->NPixel; j++) {
    
      for (i=0; i<NBuffer->BufferSize; i++) {

        /* Create a complex white noise spectrum */
        U=(gsl_rng_uniform(r)-0.5)*pi;
	G=gsl_ran_gaussian(r,sigma); 
	in[i]=cos(U)*G*Noise->WhiteRMS + sin(U)*G*Noise->WhiteRMS*I;
	      
        /* Define the noise filter */
	/* Set initial value of zeros and poles */
        Ze = 1.0 + 0.0 * I;  /* Zeros initialisation */
        Po = 1.0 + 0.0 * I;  /* Poles initialisation */
        
	/* Calculate frequency and angular frequency for spectrum */
	f=(i+1)*(*SampFreq)/(NBuffer->BufferSize);
	w=2.0*pi*f;
	
        /* Multiply all the zeros */
        for (k=0;k<Noise->Nz;k++) {
          Ze = Ze * (1.0 + Noise->Zeros[k] * w * I);
        }
	
	/* Multiply all the poles */
	for (k=0;k<Noise->Np;k++) {
	  Po = Po * (1.0 + Noise->Poles[k] * w * I);
	}
	
	/* Calculate the filter amplitude (complex) */
	H = Noise->H0 * Ze / Po; 
	
	/* Multiply the noise filter with the white noise */
	in[i]=in[i] * sqrt(1.0/(H * conj(H)));
      }
      in[0]=0.0 + 0.0*I;
      p=(fftw_plan) fftw_plan_dft_1d(NBuffer->BufferSize,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
    
      fftw_execute(p);
    
      for (i=0;i<NBuffer->BufferSize;i++) {
        NBuffer->Buffer[i][j]=creal(out[i]) * sqrt(*SampFreq/NBuffer->BufferSize) ;
      }
    
    }
    fftw_destroy_plan(p);
    fftw_free(in); 
    fftw_free(out);
    
    return *status;
}


int destroyNoiseSpectrum(NoiseSpectrum* Noise, 
			 int* const status) 
{
    
  if(Noise->Poles!=NULL){
    free(Noise->Poles);
  }
  if(Noise->Zeros!=NULL){
    free(Noise->Zeros);
  }
  if(Noise!=NULL){
    free(Noise);
  }
    
    return *status;
}

int destroyNoiseBuffer(NoiseBuffer* NBuffer, 
		       int* const status) {
    int i;
    
    if(NBuffer!=NULL){
      if(NBuffer->Buffer!=NULL){
	for (i=0;i<NBuffer->BufferSize;i++) {
	  if(NBuffer->Buffer[i]!=NULL){
	    free(NBuffer->Buffer[i]);
	  }
	}
	free(NBuffer->Buffer);
      }
      free(NBuffer);
    }

    return *status;
}
