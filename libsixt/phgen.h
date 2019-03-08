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

#ifndef PHGEN_H
#define PHGEN_H 1

#include "sixt.h"
#include "attitude.h"
#include "gendet.h"
#include "photon.h"
#include "sourcecatalog.h"


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


int phgen(Attitude* const ac,
	  SourceCatalog** const srccat,
	  const unsigned int ncat,
	  const double t0,
	  const double tend,
	  const double mjdref,
	  const double dt,
	  const float fov,
	  Photon* const ph,
	  int* const status);


#endif /* PHGEN_H */
