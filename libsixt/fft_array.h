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
