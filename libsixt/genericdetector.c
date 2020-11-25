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

#include "genericdetector.h"


int initGenericDetector(GenericDetector* gd,
			struct GenericDetectorParameters* parameters)
{
  int status = EXIT_SUCCESS;

  // Set the charge cloud dimensions.
  // Gaussian:
  gd->gcc.ccsigma =    parameters->ccsigma;
  gd->gcc.ccsize  = 3.*parameters->ccsigma;
  headas_chat(5, "Split events: assuming Gaussian Charge Cloud model:\n");
  headas_chat(5, " charge cloud size: %lf mum (3 sigma)\n",
	      gd->gcc.ccsize*1.e6);
  // Exponential model for eROSITA according to Konrad Dennerl:
  gd->ecc.parameter = 0.355;

  // Set the event thresholds:
  gd->pha_threshold = parameters->pha_threshold;
  gd->energy_threshold = parameters->energy_threshold;

  // Read the detector RMF and EBOUNDS from the specified file and
  // assign them to the Detector data structure.
  gd->rmf = loadNormalizedRMF(parameters->rmf_filename, &status);
  if(EXIT_SUCCESS!=status) return(status);

  if (strlen(parameters->arf_filename)>0) {
    gd->arf = loadARF(parameters->arf_filename, &status);
    if(EXIT_SUCCESS!=status) return(status);
  }

  return(status);
}


void cleanupGenericDetector(GenericDetector* gd)
{
  // Actually one has to release the energy allocated for the internal
  // objects in the RMF data structure. But unfortunately there exists
  // no routine with this functionality in the HEASP library.

  if (NULL==gd->rmf) return;
  freeRMF(gd->rmf);
  gd->rmf=NULL;
}


double gaussint(const double x)
{
  return(gsl_sf_erf_Q(x));
}
