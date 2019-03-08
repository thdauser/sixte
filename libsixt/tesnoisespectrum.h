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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef TES_NOISESPECTRUM_H
#define TES_NOISESPECTRUM_H 1

#include <math.h>
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

///** Noise spectrum input parameters */
//typedef struct {
//  /** White noise RMS value */
//  double WhiteRMS;
//
//  /** Normalisation of the filter function */
//  double H0;
//
//  /** Number of zeros */
//  int Nz;
//
//  /** Zeros */
//  double *Zeros;
//
//  /** Number of poles */
//  int Np;
//
//  /** Poles */
//  double *Poles;
//
//} NoiseSpectrum;


/** Noise buffer (output) */
typedef struct {
  /** Size of buffer */
  int BufferSize;

  /** Number of Pixels (to be obtained from other struct later) */
  int NPixel;

  /** Actual buffer */
  double **Buffer;
} NoiseBuffer;


/** 1/f noise generation */
typedef struct {
  /** Precision of 1/f noise generation depends on the length of
      the random-value arrays (RValues). This is the correlation 'length'. */
  int Length;

  /** Noise level */
  double Sigma;

  /** Array of random values from Gauss distribution */
  double *RValues;

  /** Sum of the random values */
  double Sumrval;

  /** Index of the array element to be changed */
  int Index;

} NoiseOoF;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////

/** Function to generate a rng for the noise buffer */
void setNoiseGSLSeed(gsl_rng **r, unsigned long int seed);

///** Function to allocate and fill noise parameter struct */
//NoiseSpectrum* newNoiseSpectrum(AdvDet *det,
//				int* const status);
NoiseBuffer* newNoiseBuffer(int* const status,
			    int *NumberOfPixels);

/** Function to initialise arrays for 1/f noise generation */
NoiseOoF* newNoiseOoF(int* const status,gsl_rng **r,double sample_freq,AdvPix* pixel);

/** Generate noise data from a noise spectrum */
int genNoiseSpectrum(AdvPix** simulated_pixels,
		     NoiseBuffer* NBuffer,
		     double *SampFreq,
		     gsl_rng **r,
		     int* const status);

void getNextOoFNoiseSumval(NoiseOoF** OFNoise,  /* */
                          gsl_rng **r,        /* Random number generator */
                          int Nactive);

/** Function to deallocate buffer memory */
int destroyNoiseBuffer(NoiseBuffer* NBuffer,
		       int* const status);
//int destroyNoiseSpectrum(NoiseSpectrum* Noise,
//			 int* const status);

int destroyNoiseOoF(NoiseOoF* OFNoise,
		       int* const status);

#endif /* TES_NOISESPECTRUM_H */
