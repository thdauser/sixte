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

fftw_complex* FFTOfArray_1d(double* Image1d, int ImageSize1, int ImageSize2);

fftw_complex* FFTOfArrayInverse_1d(double* Image1d, int ImageSize1, int ImageSize2);

fftw_complex* FFTOfArrayInverse_fftwcomplex(fftw_complex* Input, int ImageSize1, int ImageSize2);

void testFitsImagefft(fftw_complex* Image,char filename, int Size1, int Size2);

void testFitsImage1d(double* Image1d,char filename, int Size1, int Size2);

#endif
