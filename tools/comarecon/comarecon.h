#ifndef COMARECON_H
#define COMARECON_H 1


#include "sixt.h"
#include "comaevent.h"
#include "comaeventfile.h"
#include "squarepixels.h"
#include "codedmask.h"
#include "sourceimage.h"


#define TOOLSUB comarecon_main
#include "headas_main.c"


struct Parameters {
  char mask_filename[MAXMSG]; // input: coded mask reconstruction file
  char eventlist_filename[MAXMSG];
  char eventlist_template[MAXMSG];
  char image_filename[MAXMSG]; // output: reconstructed source image

  /** Detector width in [pixel]. */
  int width;
  /** Width of one detector pixel in [m]. */
  double pixelwidth;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Reads the program parameters using PIL. */
int comarecon_getpar(struct Parameters* parameters);


#endif /* COMARECON_H */
