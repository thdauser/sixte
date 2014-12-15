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


   Copyright 2007-2014 Christian Schmid, Mirjam Oertel, FAU
*/

#ifndef COMAIMGPM_H
#define COMAIMGPM_H 1


#include "sixt.h"
#include "vector.h"
#include "check_fov.h"
#include "photon.h"
#include "photonfile.h"
#include "telescope.h"
#include "codedmask.h"
#include "impact.h"
#include "impactfile.h"
#include "attitude.h"
#include "attitudefile.h"
#include "det_phi_max.h"


#define TOOLSUB comaimgPM_main
#include "headas_main.c"


struct Parameters {
  char PhotonList[MAXMSG]; //input:photon list
  char Mask[MAXMSG];       //input:mask
  char ImpactList[MAXMSG]; //output:impact list
  char Attitude[MAXMSG];   //input:attitude
  
  //Distance between the mask and detection plane ([m])
  float MaskDistance;
  //Distance between mask and collimator ([m])
  float CollDistance;
  //detector-size (width, depth) ([m])
  float x_det, y_det;
  //width of one detector pixel ([m])
  float det_pixelwidth;
  //Width of one detector-element that is surrounded by the collimator ([m])
  float det_width;
  //Width of collimator wall ([m])
  float wall;

  //detector pointing direction
  double RA, DEC;

  //time-offset ([s])
  double Timezero;

  //Exposure time ([s])
  double Exposure;
};


// Function declarations.

/** Reads the program parameters using PIL. */
int comaimgPM_getpar(struct Parameters* par);


#endif /* COMAIMGPM_H */
