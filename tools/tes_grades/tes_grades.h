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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef TES_GRADES_H
#define TES_GRADES_H 1

#include "sixt.h"
#include "eventfile.h"


#define TOOLSUB tes_grades_main
#include "headas_main.c"


struct Parameters{
  char EventList[MAXFILENAME];

  /** Characteristic time unit of the TES microcalorimeter. */
  double TimeUnit;
  int PreTrigger;
  int PostTrigger;
  double PileupTime;

  char clobber;
};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


int tes_grades_getpar(struct Parameters* par);


#endif /* TES_GRADES_H */
