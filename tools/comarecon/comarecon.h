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
#include "find_position.h"
#include "vector.h"
#include "maskshadow.h"


#define TOOLSUB comarecon_main
#include "headas_main.c"


struct Parameters {
  char Mask[MAXMSG]; // input: coded mask reconstruction file
  char EventList[MAXMSG];
  char EventListTemplate[MAXMSG];
  char Image[MAXMSG]; // output: reconstructed source image
  char PositionList[MAXMSG]; // output: table of found sources

  //detector pointing direction
  double RA, DEC;

  /** Detector width in [pixel]. */
  int width;
  /** Width of one detector pixel in [m]. */
  double pixelwidth;
  /** Width of re-pixeled detector pixels in [m]. */
  double RePixSize;  //0.0 means detector shouldn't be repixeled, else the value is given

  /** Distance between the coded mask and the detector plane ([m]). */
  double MaskDistance; 

  /**length of DCU, gap between 2 DCU's and gap between two DCA's [m]. */
  //only works for 2x2 DCU's separated by DCU_gap, followed by DCA_gap
  double DCU_length;
  double DCU_gap;
  double DCA_gap;
  
  /**threshold for sources, factor to mulpilpy sigma with. */
  double Sigma;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int comarecon_getpar(struct Parameters* parameters);


#endif /* COMARECON_H */
