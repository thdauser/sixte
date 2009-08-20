#ifndef PSF_H
#define PSF_H 1

#include "sixt.h"

#include "telescope.h"
#include "photon.h"
#include "point.h"
#include "vignetting.h"


/** Stores the PSF data for one particular off-axis angle 
 * and one particular energy. */
typedef struct {
  double** data;   /**< pointer to PSF data array [x][y]. */

  double angle;    /**< off-axis angle of this particular PSF [rad]. */
  double energy;   /**< energy of this particular PSF [keV]. */

  int naxis1, naxis2;    /**< Width of the image [pixel]. */
  double cdelt1, cdelt2; /**< Width of one pixel [rad]. */
  double crpix1, crpix2; /**< [pixel] */
  double crval1, crval2; /**< [rad] */
} PSF_Item;


/** Storage for the several PSFs available for a mirror system. */
typedef struct {
  /** Number of PSF_Items in this store (#(offaxis-angles)*#(energies)). */
  int N_elements;

  /** Array of PSF_Items for the different discrete off-axis angles and energies. */
  PSF_Item *item;
} PSF;


///////////////////////////////////////////////////////////////////////////////
// Function declarations
///////////////////////////////////////////////////////////////////////////////


// Reads the specified PSF FITS file (containing images) and
// stores these data in the PSF storage. 
PSF* get_psf(const char*, int* status);


/** Calculates the position on the detector, where a photon at given sky position 
 * with specified energy hits the detector according to the PSF data. 
 * The exact position is determined with a random number generator 
 * (randomization over one PSF pixel).
 * Return value is '1', if the photon hits the detector. If it does not fall onto the
 * detector, the function returns '0'.  The output detector position is stored 
 * in [m] in the first 2 parameters of the function. */
int get_psf_pos(struct Point2d*, Photon, struct Telescope, Vignetting*, PSF*);


/** Releases the memory of the PSF storage. */
void free_psf(PSF*);


// Save the data contained in the PSF storage to images in a FITS file
// according to the OGIP standards.
int save_psf_image(PSF*, const char*, int* status);


#endif /* PSF_H */

