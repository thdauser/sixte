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


   Copyright 2007-2014 Christian Schmid, FAU
*/

#ifndef XMS_PIXTEMP_H
#define XMS_PIXTEMP_H 1

#include "sixt.h"
#include "event.h"
#include "eventfile.h"
#include "rmf.h"

#define TOOLSUB xms_pixtemp_main
#include "headas_main.c"


struct Parameters{
  /** Filename of the XMS event file. */
  char EventList[MAXMSG];

  /** Filename of the output file. */
  char OutputFile[MAXMSG];

  /** Filename of the detector response file containing the EBOUDNS
      table. */
  char RSP[MAXMSG];

  /** X- and y- coordinate of the pixel to be investigated. */
  int pixx, pixy;

};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


int xms_pixtemp_getpar(struct Parameters*);


#endif /* XMS_PIXTEMP_H */
