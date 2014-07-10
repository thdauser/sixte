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

#ifndef TESSTREAM_H
#define TESSTREAM_H 1

#include "sixt.h"
#include "gti.h"
#include "pixelimpact.h"
#include "advdet.h"
#include "pixelimpactfile.h"
#include "tesproftemplates.h"
#include "tesnoisespectrum.h"
#include "tesdatastream.h"

#define TOOLSUB tesstream_main
#include "headas_main.c"

struct Parameters {
  char PixImpList[MAXFILENAME];
  char XMLFile[MAXFILENAME];
  char streamname[MAXFILENAME];
  
  char activePixels[9];
  int Nactive;
  int Npix;
  int nlo;
  int nhi;
  
  double tstart;
  double tstop;
  
  char clobber;
  char history;
  
  unsigned long int seed;
};

int getpar(struct Parameters* const par);

#endif /* TESSTREAM_H */