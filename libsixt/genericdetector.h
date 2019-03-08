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

#ifndef GENERICDETECTOR_H
#define GENERICDETECTOR_H 1

#include "sixt.h"
#include "gaussianchargecloud.h"
#include "exponentialchargecloud.h"
#include "rmf.h"
#include "arf.h"


/** The GNU Scientific Library Errorfunction is used to calculate
    charge distribution of split events (assuming a Gaussian shape for
    the carge cloud). */
#include <gsl/gsl_sf_erf.h>

/** This value represents an invalid pixel. */
#define INVALID_PIXEL (-1)


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic x-ray detector. Contains common data/specifications that
    are defined for all different kinds of detectors.  The
    GenericDetector data structure is usually included as element of
    more specific detector models like, e.g. WFIDetector or
    FramestoreDetector.*/
typedef struct {

  /** Properties of Gaussian shaped charge clouds. */
  GaussianChargeCloud gcc;

  /** Properties of exponential charge clouds. */
  ExponentialChargeCloud ecc;

  /** Lower detector PHA threshold [PHA channels]. Events with a
      lower PHA value are dismissed. If the PHA threshold is set to
      "-1" the energy threshold is taken into account. */
  long pha_threshold;
  /** Lower detector energy threshold [keV]. Events with a lower
      energy / pixel charge are dismissed. This value is only
      regarded if the pha_threshold is set to "-1". */
  float energy_threshold;

  /** Detector response matrix. The RSP file that is originally loaded
      may also contain ARF contributions. But as they already have
      been taken into account in the generation of the input spectra
      for the X-ray sources, the ARF contributions have to be removed
      by normalizing the RSP matrix to obtain an RMF. */
  struct RMF* rmf;

  /** Detector ARF. */
  struct ARF* arf;

} GenericDetector;


/** Data required to initialize the GenericDetector data structure.
    The meaning of the inidividual parameters is the same as in the
    GenericDetector data structure. */
struct GenericDetectorParameters {
  double ccsigma;

  long pha_threshold;
  float energy_threshold;

  char* rmf_filename;
  char* arf_filename;
};


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Set up the initial configuraton for the GenericDetector data
    structure.  This routines sets the values of the GenericDetector
    elements according to the given parameters. */
int initGenericDetector(GenericDetector*, struct GenericDetectorParameters*);

/** Clean up the GenericDetector data structure. Release the allocated
    memory. */
void cleanupGenericDetector(GenericDetector* gd);

/** Calculates the Gaussian integral using the GSL complementary error
    function. */
double gaussint(const double x);


#endif /* GENERICDETECTOR_H */
