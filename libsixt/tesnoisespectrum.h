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

#ifndef TES_NOISESPECTRUM_H 
#define TES_NOISESPECTRUM_H 1

#include <complex.h>
#include "sixt.h"
#include "advdet.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fftw3.h>

#define NOISEBUFFERSIZE 65536


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////

/** Noise spectrum input parameters */
typedef struct {
  /** White noise RMS value */
  double WhiteRMS;
  
  /** Normalisation of the filter function */
  double H0;
  
  /** Number of zeros */
  int Nz;
  
  /** Zeros */
  double *Zeros;
  
  /** Number of poles */
  int Np;
  
  /** Poles */
  double *Poles;
  
} NoiseSpectrum;


/** Noise buffer (output) */
typedef struct {
  /** Size of buffer */
  int BufferSize;
  
  /** Number of Pixels (to be obtained from other struct later) */
  int NPixel;
  
  /** Actual buffer */
  double **Buffer;
} NoiseBuffer;



/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////

/** Function to allocate and fill noise parameter struct */
NoiseSpectrum* newNoiseSpectrum(AdvDet *det, 
				int* const status);
NoiseBuffer* newNoiseBuffer(int* const status, 
			    int *NumberOfPixels);


/** Generate noise data from a noise spectrum */
int genNoiseSpectrum(NoiseSpectrum* Noise, 
		     NoiseBuffer* NBuffer, 
		     double *SampFreq, 
		     int* const status);


/** Function to deallocate buffer memory */
int destroyNoiseBuffer(NoiseBuffer* NBuffer, 
		       int* const status);
int destroyNoiseSpectrum(NoiseSpectrum* Noise, 
			 int* const status);


#endif /* TES_NOISESPECTRUM_H */
