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

#ifndef PHIMG_H
#define PHIMG_H 1

#include "sixt.h"
#include "attitude.h"
#include "check_fov.h"
#include "gentel.h"
#include "photon.h"
#include "impact.h"


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


int phimg(const GenTel* const tel,
	  Attitude* const ac,
	  Photon* const ph,
	  Impact* const imp,
	  int* const status);


#endif /* PHIMG_H */
