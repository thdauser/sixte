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

#include "phgen.h"


int phgen(Attitude* const ac,
	  SourceCatalog** const srccat,
	  const unsigned int ncat,
	  const double t0,
	  const double tend,
	  const double mjdref,
	  const double dt,
	  const float fov,
	  Photon* const ph,
	  int* const status)
{
  // Photon list buffer.
  static LinkedPhoListElement* pholist=NULL;
  // Counter for the photon IDs.
  static long long ph_id=0;

  // Current time.
  static double time=0.;
  if (time<t0) {
    time=t0;
  }

  // If the photon list is empty generate new photons from the
  // given source catalog.
  while((NULL==pholist)&&(time<tend)) {
    // Determine the telescope pointing at the current point of time.
    Vector pointing=getTelescopeNz(ac, time, status);
    CHECK_STATUS_BREAK(*status);

    // Display the program progress status.
    double ra, dec;
    calculate_ra_dec(pointing, &ra, &dec);

    // Generate new photons for all specified catalogs.
    double t1=MIN(time+dt, tend);
    unsigned int ii;
    for (ii=0; ii<ncat; ii++) {
      if (NULL==srccat[ii]) continue;

      // Get photons for all sources in the catalog.
      LinkedPhoListElement* newlist=
	genFoVXRayPhotons(srccat[ii], &pointing, fov,
			  time, t1, mjdref, status);
      CHECK_STATUS_BREAK(*status);

      // Merge the photon lists.
      pholist=mergeLinkedPhoLists(pholist, newlist);
    }

    // Increase the time.
    time+=dt;
  }

  // If there is no photon in the buffer.
  if (NULL==pholist) return(0);

  // Take the first photon from the list and return it.
  copyPhoton(ph, &pholist->photon);

  // Delete the processed element.
  LinkedPhoListElement* next=pholist->next;
  free(pholist);
  pholist=next;

  // Set the photon ID.
  ph->ph_id=++ph_id;

  return(1);
}
