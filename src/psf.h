#ifndef PSF_H
#define PSF_H 1


/** 
 * This file contains all definitions required for PSF calculations.
 */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>

#include "fitsio.h"
#include "headas.h"
#include "headas_error.h"


#define PSF_NFIELDS 5 // number of fields in the table (= number of columns)



//////////////////////////////////////////////////////////////////////////////


// This structure is used to store the PSF data for one particular 
// off-axis angle and one particular energy.
typedef struct {
  double** data;   // pointer to PSF data array [x][y]

  double angle;    // off-axis angle of this particular PSF [rad]
  double energy;   // energy of this particular PSF [keV]

  // For each energy the PSFs are rescaled in such a way that the on-axis PSF 
  // is normalized to 1.  The scaling factor is stored here and has to be used 
  // for rescaling the source spectra appropriately.
  double scaling_factor;
} PSF_Item;



// Storage for the several PSF data units that belong to one mirror system.
typedef struct {
  // Number of PSF arrays in this store (#(offaxis-angles)*#(energies)):
  int N_elements;
  // Width of the PSF in pixels:
  int width;     
  // Width of a single PSF pixel in [m]:
  double pixelwidth; 

  // PSF data for the individual off-axis angles and energies
  PSF_Item *item;
} PSF;



#include "detectors.h"
#include "vector.h"
#include "random.h"
#include "telescope.h"
#include "photon.h"
#include "point.h"



///////////////////////////////////////////////////////////////////////////////
// Function declarations
///////////////////////////////////////////////////////////////////////////////


// Reads the specified PSF FITS file (containing images) and
// stores these data in the PSF storage. 
PSF* get_psf(const char*, int* status);


// Calculates the position on the detector, where a photon at given sky position 
// with specified energy hits the detector according to the PSF data and a random 
// number generator (randomization over one PSF pixel).
// Return value is '1', if the photon hits the detector. If it does not fall onto the
// detector, the function returns '0'.  The output detector position is stored 
// in [mu m] in the first 2 parameters of the function.
int get_psf_pos(struct Point2d*, struct Photon, struct Telescope, PSF*);


// Releases the memory for the PSF storage.
void free_psf(PSF*);


// Save the data contained in the PSF storage to images in a FITS file
// according to the OGIP standards.
int save_psf_image(PSF*, const char *, int *status);


#endif /* PSF_H */

