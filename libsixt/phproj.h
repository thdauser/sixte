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

#ifndef PHPROJ_H
#define PHPROJ_H 1

#include "sixt.h"
#include "attitude.h"
#include "event.h"
#include "eventfile.h"
#include "geninst.h"
#include "point.h"
#include "vector.h"


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


void phproj(GenInst* const inst,
	    Attitude* const ac,
	    EventFile* const plf,
	    const double t0,
	    const double exposure,
	    int* const status);


#endif /* PHPROJ_H */
