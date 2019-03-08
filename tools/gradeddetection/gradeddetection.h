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


   Copyright 2014 Philippe Peille, IRAP
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef GRADEDDETECTION_H
#define GRADEDDETECTION_H 1

#include "sixt.h"
#include "parinput.h"
#include "advdet.h"
#include "crosstalk.h"
#include "pixelimpactfile.h"
#include "tesinitialization.h"
#include "grading.h"


#define TOOLSUB gradeddetection_main
#include "headas_main.c"

struct Parameters {
  char* PixImpList;
  char* EvtFile;
  char* AdvXml;

  double tstart;
  double tstop;

  int doCrosstalk;
  int saveCrosstalk;

  int clobber;
  int history;

  int seed;
};

void gradeddetection_getpar(struct Parameters* const par,int* const status);

#endif /* GRADEDDETECTION_H */
