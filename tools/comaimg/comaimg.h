#ifndef COMAIMG_H
#define COMAIMG_H 1


#include "sixt.h"
#include "vector.h"
#include "check_fov.h"
#include "photon.h"
#include "photonlistfile.h"
#include "telescope.h"
#include "codedmask.h"
#include "impact.h"
#include "impactlistfile.h"
#include "attitudecatalog2.h"
#include "attitude.h"

#include "det_phi_max.h"


#define TOOLSUB comaimg_main
#include "headas_main.c"


struct Parameters {
  char PhotonList[MAXMSG]; //input:photon list
  char Mask[MAXMSG];       //input:mask
  char ImpactList[MAXMSG]; //output:impact list
  char Attitude[MAXMSG];   //input:attitude
  
  //Distance between the mask and detection plane ([m])
  float MaskDistance;
  //mask-size (width, depth) ([m])
  float x_mask, y_mask;
  //detector-size (width, depth) ([m])
  float x_det, y_det;
  //width of one detector pixel ([m])
  float det_pixelwidth;

  //detector pointing direction
  double RA, DEC;

  //time-offset ([s])
  double Timezero;

  //Exposure time ([s])
  double Exposure;
};


// Function declarations.

/** Reads the program parameters using PIL. */
int comaimg_getpar(struct Parameters* par);


#endif /* COMAIMG_H */
