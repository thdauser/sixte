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

void setNoiseGSLSeed(gsl_rng **r, unsigned long int seed){
  
  const gsl_rng_type * T;
  T = gsl_rng_default;
  *r = gsl_rng_alloc(T);
  gsl_rng_set(*r, seed);
}


//NoiseSpectrum* newNoiseSpectrum(AdvDet *det,
//				int* const status) {
//
//    NoiseSpectrum* Noise=malloc(sizeof(NoiseSpectrum));
//    if(Noise==NULL){
//      *status=EXIT_FAILURE;
//      SIXT_ERROR("memory allocation for noise spectrum failed");
//      CHECK_STATUS_RET(*status, Noise);
//    }
//    int ii;
//
//    /* Set noise parameters */
//    Noise->WhiteRMS=det->TESNoise->WhiteRMS;  /* White Noise normalisation */
//    Noise->H0=det->TESNoise->H0;              /* Noise filter normalisation */
//    Noise->Nz=det->TESNoise->Nz;              /* Number of zeros in the noise spectrum */
//    Noise->Np=det->TESNoise->Np;              /* Number of poles in the noise spectrum */
//
//    Noise->Poles=NULL;
//    Noise->Zeros=NULL;
//
//    if(Noise->Np>0){
//      Noise->Poles=(double*)malloc(Noise->Np*sizeof(double));
//      if(Noise->Poles==NULL){
//	*status=EXIT_FAILURE;
//	SIXT_ERROR("memory allocation for noise spectrum failed");
//	CHECK_STATUS_RET(*status, Noise);
//      }
//      for(ii=0; ii<Noise->Np; ii++){
//	Noise->Poles[ii]=det->TESNoise->Poles[ii];
//      }
//    }
//    if(Noise->Nz>0){
//      Noise->Zeros=(double*)malloc(Noise->Nz*sizeof(double));
//      if(Noise->Zeros==NULL){
//	*status=EXIT_FAILURE;
//	SIXT_ERROR("memory allocation for noise spectrum failed");
//	CHECK_STATUS_RET(*status, Noise);
//      }
//      for(ii=0; ii<Noise->Nz; ii++){
//	Noise->Zeros[ii]=det->TESNoise->Zeros[ii];
//      }
//    }
//    return Noise;
//}


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

NoiseOoF* newNoiseOoF(int* const status,gsl_rng **r,double sample_freq,AdvPix* pixel) {
    int j;
    double OffSet, diff;
    
    /* Allocate memory for 1/f noise arrays */
    NoiseOoF* OFNoise=(NoiseOoF*) malloc(sizeof(NoiseOoF));
    if(OFNoise==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for 1/f noise generation failed");
      CHECK_STATUS_RET(*status, OFNoise);
    } 
    
    OffSet=(double) pixel->ADCOffset;
    
    /* Length of the 1/f noise arrays */
    OFNoise->Length=(int) (sample_freq * 0.3678)/pixel->TESNoise->OoFKnee;
    /* Noise level */
    OFNoise->Sigma=pixel->TESNoise->OoFRMS;
    
    OFNoise->RValues=(double*)malloc(OFNoise->Length*sizeof(double));
    
    if(OFNoise->RValues==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for OFNoise RValues failed");
      CHECK_STATUS_RET(*status, OFNoise);
    }
    
    OFNoise->Sumrval=0.;
    for(j=0;j<OFNoise->Length;j++) {
        /* Fill array with Gauss-distributed random values and sum the array */
	OFNoise->RValues[j]=gsl_ran_gaussian(*r,OFNoise->Sigma);
	OFNoise->Sumrval=OFNoise->Sumrval+OFNoise->RValues[j];
    }	
    /* Make sure the data is not too far away from the baseline */
    diff=abs(OFNoise->Sumrval);
    while (diff >= OffSet/5) {
	OFNoise->Sumrval=0.;
	for(j=0;j<OFNoise->Length;j++) {
	  OFNoise->RValues[j]=gsl_ran_gaussian(*r,OFNoise->Sigma);
	  OFNoise->Sumrval=OFNoise->Sumrval+OFNoise->RValues[j];
	}  
	diff=abs(OFNoise->Sumrval);
    }
    
    OFNoise->Index=0;
    
    
    return OFNoise;
}


int genNoiseSpectrum(AdvPix** simulated_pixels,
		     NoiseBuffer* NBuffer, 
		     double *SampFreq, 
		     gsl_rng **r,
                     int* const status) 
{
    double Gx, Gy, sigma, df;
    const double pi=3.14159265359;
    int i,j, k;
    
    fftw_complex *in;
    double *out;
    fftw_plan p;
    fftw_complex Ze, Po, H; /* Products of Zeros & Poles */
    double w, f;               /* Omega and frequency*/
    
    sigma=1.;
    
    in=(fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (NBuffer->BufferSize/2+1));
    if(in==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for fftw_complex failed");
      CHECK_STATUS_RET(*status, *status);
    }
    out=(double *) fftw_malloc(sizeof(double) * NBuffer->BufferSize);
    if(out==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for fftw_complex failed");
      CHECK_STATUS_RET(*status, *status);
    }
    
    
    /* Calculate size of frequency bin */
    df=*SampFreq/(NBuffer->BufferSize);

    for (j=0; j<NBuffer->NPixel; j++) {
    
      for (i=1; i<=NBuffer->BufferSize/2; i++) {
        
        /* Define the noise filter */
	/* Set initial value of zeros and poles */
        Ze = 1.0 + 0.0 * I;  /* Zeros initialisation */
        Po = 1.0 + 0.0 * I;  /* Poles initialisation */
        
	/* Calculate frequency and angular frequency for spectrum */
	f=i*df;
	w=2.0*pi*f;

	/* Multiply all the zeros */
        for (k=0;k<simulated_pixels[j]->TESNoise->Nz;k++) {
          Ze = Ze * (1.0 + simulated_pixels[j]->TESNoise->Zeros[k] * w * I);
        }
	
	/* Multiply all the poles */
	for (k=0;k<simulated_pixels[j]->TESNoise->Np;k++) {
	  Po = Po * (1.0 + simulated_pixels[j]->TESNoise->Poles[k] * w * I);
	}
	
	/* Calculate the filter amplitude (complex) */
	H = simulated_pixels[j]->TESNoise->H0 * Ze / Po;
	
	if (i==NBuffer->BufferSize/2) { // At Nyquist freq, the FT is purely real-> draw only one gaussian variable
	  /* Create a complex white noise spectrum */
	  Gx=gsl_ran_gaussian(*r,sigma); 
	  in[i]=Gx + 0.0*I;
	  
	  /* Multiply the noise filter with the white noise */
	  in[i]=in[i] * abs(H) * simulated_pixels[j]->TESNoise->WhiteRMS * sqrt(df) / sqrt(2.);
	} else {
	  /* Create a complex white noise spectrum */
	  Gx=gsl_ran_gaussian(*r,sigma); 
	  Gy=gsl_ran_gaussian(*r,sigma); 
	  in[i]=Gx + Gy*I;

	  /* Multiply the noise filter with the white noise */
	  in[i]=in[i] * H * simulated_pixels[j]->TESNoise->WhiteRMS *sqrt(df) / sqrt(2.);
	}
      }
      in[0]=0.0 + 0.0*I;
      p=(fftw_plan) fftw_plan_dft_c2r_1d(NBuffer->BufferSize,in,out,FFTW_ESTIMATE);
    
      fftw_execute(p);
    
      for (i=0;i<NBuffer->BufferSize;i++) {
        NBuffer->Buffer[i][j]=out[i] / sqrt(2*NBuffer->BufferSize) ;
      }
    
    }
    fftw_destroy_plan(p);
    fftw_free(in); 
    fftw_free(out);
    
    return *status;
}

void getNextOoFNoiseSumval(NoiseOoF** OFNoise,  /* */
                          gsl_rng **r,        /* Random number generator */
                          int Nactive)
{
	for(int pixNumber=0;pixNumber<Nactive;pixNumber++){
		if (OFNoise[pixNumber]!=NULL){
			int i, c;
			double src;

			i=OFNoise[pixNumber]->Index;

			i++;
			OFNoise[pixNumber]->Index=i;

			if (i==OFNoise[pixNumber]->Length) {
				OFNoise[pixNumber]->Index=0;
				continue;
			}

			c=0;
			while( (i & 1) == 0 ) {
				i = i >> 1;
				c++;
			}

			src=gsl_ran_gaussian(*r,OFNoise[pixNumber]->Sigma);
			OFNoise[pixNumber]->Sumrval=OFNoise[pixNumber]->Sumrval + src - OFNoise[pixNumber]->RValues[c];
			OFNoise[pixNumber]->RValues[c]=src;
		}
	}
}		      


//int destroyNoiseSpectrum(NoiseSpectrum* Noise,
//			 int* const status)
//{
//
//  if(Noise->Poles!=NULL){
//    free(Noise->Poles);
//  }
//  if(Noise->Zeros!=NULL){
//    free(Noise->Zeros);
//  }
//  if(Noise!=NULL){
//    free(Noise);
//  }
//
//    return *status;
//}

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

int destroyNoiseOoF(NoiseOoF* OFNoise, 
		    int* const status) {

    if(OFNoise!=NULL){
      if(OFNoise->RValues!=NULL){
	free(OFNoise->RValues);
      }
      free(OFNoise);
    }

    return *status;
}

