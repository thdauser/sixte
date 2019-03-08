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

#ifndef VIGNETTING_H
#define VIGNETTING_H 1

#include "sixt.h"


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** Data structure containing the mirror vignetting function. */
typedef struct {
  /** Number of energy bins. */
  int nenergies;
  /** Number of off-axis angles. */
  int ntheta;
  /** Number of azimuth angles. */
  int nphi;

  /** Mean energy of bin in [keV] (0.5*(energ_lo+energ_hi). */
  float* energy;
  /** Off-axis angle in [rad]. */
  float* theta;
  /** Azimuthal angle in [rad]. */
  float* phi;
  /** Vignetting data. Array[energy, theta, phi] */
  float*** vignet;

  /** Minimum available energy [keV]. */
  float Emin;
  /** Maximum available energy [keV]. */
  float Emax;

} Vignetting;


///////////////////////////////////////////////////////////////////////////////
// Function declarations.
///////////////////////////////////////////////////////////////////////////////


/** Constructor of the Vignetting data structure. Loads the
    vignetting function from a given FITS file. The format of the
    FITS file is defined by OGIP Memo CAL/GEN/92-021. */
Vignetting* newVignetting(const char* const filename, int* const status);

/** Destructor for Vignetting data structure. */
void destroyVignetting(Vignetting** const vi);

/** Determine the Vignetting factor for given photon energy, off-axis
    angle, and azimuth angle. The energy has to be given in [keV], the
    angles in [rad]. If the pointer to the Vignetting data structure
    is NULL, a default value of 1. will be returned. Note that the
    azimuthal dependence is neglected in the current
    implementation. */
float get_Vignetting_Factor(const Vignetting* const vi,
			    const float energy,
			    const float theta,
			    const float phi);


#endif /* VIGNETTING_H */
