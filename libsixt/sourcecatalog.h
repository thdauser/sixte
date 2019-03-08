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

#ifndef SOURCECATALOG_H
#define SOURCECATALOG_H 1

#include "sixt.h"
#include "check_fov.h"
#include "gendet.h"
#include "kdtreeelement.h"
#include "linkedpholist.h"
#include "simput.h"
#include "source.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Catalog of different X-ray sources. */
typedef struct {
  /** KDTree containing the Source objects for all point-like
      sources. */
  KDTreeElement* tree;

  /** Array containing Source objects for all extended sources. */
  Source* extsources;

  /** Number of entries in the linear array of extended sources. */
  long nextsources;

  /** SIMPUT source catalog containing all relevant data. */
  SimputCtlg* simput;

} SourceCatalog;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. */
SourceCatalog* newSourceCatalog(int* const status);

/** Destructor. */
void freeSourceCatalog(SourceCatalog** const cat, int* const status);

/** Load a SIMPUT source catalog from a FITS file. */
SourceCatalog* loadSourceCatalog(const char* const filename,
				 struct ARF* const arf,
				 int* const status);

/** Create photons for all sources in the catalog for the specified
    time interval. Only sources within the FoV (diameter given in
    [rad]) around the telescope pointing direction are taken into
    account. */
LinkedPhoListElement* genFoVXRayPhotons(SourceCatalog* const cat,
					const Vector* const pointing,
					const float fov,
					const double t0, const double t1,
					const double mjdref,
					int* const status);


#endif /* SOURCECATALOG_H */
