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


   Copyright 2007 - 2018: Christian Schmid, Mirjam Oertel, FAU.
   Manuel Castro, National Institute for Space Research (INPE),
		 Brazil; under grant #2017/00968-6,
		 São Paulo Research Foundation (FAPESP).
   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/


#ifndef RUNMASK_H
#define RUNMASK_H 1

#include "sixt.h"

#include "attitude.h"
#include "geninst.h"
#include "phgen.h"
#include "photonfile.h"
#include "simput.h"
#include "vector.h"
#include "sourcecatalog.h"
#include "impact.h"
#include "impactfile.h"
#include "comadetector.h"
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
#include "maskshadow.h"
#include "testimg.h"
#include "repix.h"
#include "check_fov.h"
#include "photon.h"
#include "telescope.h"
#include "attitudefile.h"
#include "det_phi_max.h"
#include "masksystem.h"


#define TOOLSUB runmask_main
#include "headas_main.c"


struct Parameters {
  char PhotonList[MAXFILENAME];
  char Mission[MAXMSG];
  char Instrument[MAXMSG];
  char Mode[MAXMSG];
  char XMLFile[MAXFILENAME];
  char Attitude[MAXFILENAME];
  char ImpactList[MAXFILENAME];
  char EventList[MAXMSG];
  char EventListTemplate[MAXMSG];
  char Image[MAXMSG];
  char PositionList[MAXMSG];
  char AdvXMLFile[MAXFILENAME]; //Advanced setup for coded mask systems

    //detector pointing direction
  double RA, DEC;

  char Simput[MAXFILENAME];

  double MJDREF;
  double TSTART;
  double dt; //---> ???

//time-offset ([s])
  double Timezero;

//Exposure time ([s])
  double Exposure;

/**threshold for sources, factor to mulpilpy sigma with. */
  double Sigma;

  int Seed;

  char clobber;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


int runmask_getpar(struct Parameters* parameters);
int photogen(struct Parameters* parameters, GenInst* inst);
int comaimg(struct Parameters* parameters, MaskSystem* mask_setup);
int comadet(struct Parameters* parameters, MaskSystem* mask_setup);
int comarecon(struct Parameters* parameters, MaskSystem* mask_setup);

#endif /* */
