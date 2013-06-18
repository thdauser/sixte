#ifndef COMARECON_H
#define COMARECON_H 1


#include "sixt.h"
#include "comaevent.h"
#include "comaeventfile.h"
#include "squarepixels.h"
#include "codedmask.h"
#include "sourceimage.h"
#include "reconstruction.h"
#include "eventarray.h"
#include "fft_array.h"
#include "fftw3.h"
#include "balancing.h"


#define TOOLSUB comarecon_main
#include "headas_main.c"


struct Parameters {
  char Mask[MAXMSG]; // input: coded mask reconstruction file
  char EventList[MAXMSG];
  char EventListTemplate[MAXMSG];
  char Image[MAXMSG]; // output: reconstructed source image
  char ReconArray[MAXMSG]; //output: reconstruction array

  int ReconType; // reconstruction type (1: direct deconvolution, 2: FFT)

  /** Detector width in [pixel]. */
  int width;
  /** Width of one detector pixel in [m]. */
  double pixelwidth;

  /** Distance between the coded mask and the detector plane ([m]). */
  double MaskDistance;  
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int comarecon_getpar(struct Parameters* parameters);


#endif /* COMARECON2_H */
