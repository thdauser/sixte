#ifndef VIGNETTING_H
#define VIGNETTING_H 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <assert.h>
#include "fitsio.h"

#include "sixt.h"
#include "headas.h"
#include "headas_error.h"


/** Data structure containing the mirror vignetting function. */
typedef struct {
  int nenergies; /**< Number of energy bins. */
  int ntheta; /**< Number of off-axis angles. */
  int nphi; /**< Number of azimuth angles. */

  float* energ_lo; /**< Minimum energy of bin in [keV]. */
  float* energ_hi; /**< Maximum energy of bin in [keV]. */
  float* theta; /**< Off-axis angle in [rad]. */
  float* phi; /**< Azimuth angle in [rad]. */
  float*** vignet; /**< Vignetting data. Array[energy, theta, phi] */

  float Emin; /**< Minimum available energy. */
  float Emax; /**< Maximum available energy. */
} Vignetting;



/** Constructor of the Vignetting data structure. 
 * Loads the vignetting function from a given FITS file. 
 * The format of the FITS file is defined by 
 * OGIP Memo CAL/GEN/92-021. */
Vignetting* get_Vignetting(char* filename, int* status);


/** Destructor for Vignetting data structure. */
void free_Vignetting(Vignetting* vi);


/** Determine the Vignetting factor for given photon energy, off-axis angle, 
 * and azimuth angle. 
 * The energy has to be given in  [keV], the angles in [rad]. */
/* TODO: So far the azimuth angle is neglected! */
float get_Vignetting_Factor(Vignetting* vi, float energy, float theta, float phi);


#endif // #ifndef VIGNETTING_H
