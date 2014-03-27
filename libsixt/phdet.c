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

#include "phdet.h"


void phdetGenDet(GenDet* const det,
		 Impact* const impact,
		 const double tend,
		 int* const status)
{
  double operation_time;
  if (NULL==impact) {
    // If no impact has been given as parameter, finalize the GenDet. 
    // Perform the time-triggered operations without adding any new 
    // signal charges.
    operation_time=tend;
  } else {
    operation_time=impact->time;
  }

  // Call the detector operating clock routine.
  operateGenDetClock(det, operation_time, status);
  CHECK_STATUS_VOID(*status);


  // Total number of detected photons. Only the number of
  // photons absorbed by valid pixels inside the detector is
  // counted. Split events created by one photon are counted only
  // once.
  static unsigned long n_detected_photons=0;

  // Check if an impact has been given as a parameter.
  if (NULL!=impact) {
    // Add the impact to the detector array. If it is absorbed
    // by at least one valid pixel, increase the counter for
    // the number of detected photons.
    if (addGenDetPhotonImpact(det, impact, status) > 0) {
      n_detected_photons++;
    }
    CHECK_STATUS_VOID(*status);
  }
}
