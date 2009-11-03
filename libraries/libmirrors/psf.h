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

  /** Array of PSF_Items for the different discrete off-axis angles
      and energies. */
  PSF_Item *item; // Obsolete

  /** Array of PSF_Items for different photon energies, off-axis angles, and
      azimuthal angles. */
  PSF_Item*** data;

  /** Number of different energies PSF images are available for. */
  int nenergies;
  /** Different energies PSF images are available for ([kev]). */
  double* energies;
  /** Number of different off-axis angles PSF images are available for. */
  int nthetas;
  /** Different off-axis angles PSF images are available for ([rad]). */
  double* thetas;
  /** Number of different azimuthal angles PSF images are available for. */
  int nphis;
  /** Different azimuthal angles PSF images are available for ([rad]). */
  double* phis;

} PSF;


///////////////////////////////////////////////////////////////////////////////
// Function declarations
///////////////////////////////////////////////////////////////////////////////


/** Constructor for the PSF data structure. Reads PSF data from a FITS
    file with one or several image extensions. The file format is
    given by OGIP Calibration Memo CAL/GEN/92-027. */
PSF* get_psf(const char* filename, int* status);

/** Calculates the position on the detector, where a photon at given
    sky position with specified energy hits the detector according to
    the PSF data.  * The exact position is determined with a random
    number generator * (randomization over one PSF pixel).  * Return
    value is '1', if the photon hits the detector. If it does not fall
    onto the * detector, the function returns '0'.  The output
    detector position is stored * in [m] in the first 2 parameters of
    the function. */
int get_psf_pos(/** Output: coordinates of the photon on the detector ([m]). */
		struct Point2d* position, 
		/** Incident photon. */
		Photon, 
		/** Telescope information (focal length, pointing directions. */
		struct Telescope telescope, 
		/** Vignetting function. */
		Vignetting* vignetting, 
		/** PSF with data for different photon energies and off-axis angles. */
		PSF* psf);

/** Releases the memory of the PSF storage. */
void free_psf(PSF* psf);

// Save the data contained in the PSF storage to images in a FITS file
// according to the OGIP standards.
int save_psf_image(PSF* psf, const char* filename, int* status);


#endif /* PSF_H */

