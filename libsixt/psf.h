#ifndef PSF_H
#define PSF_H 1

#include "sixt.h"

#include "telescope.h"
#include "photon.h"
#include "point.h"
#include "vector.h"
#include "vignetting.h"


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** Stores the PSF data for one particular off-axis angle and one
    particular energy. */
typedef struct {
  /** Pointer to the PSF data array [x][y]. */
  double** data;   

  /** Width of the image [pixel]. */
  int naxis1, naxis2;   
  /** Width of one pixel [m]. */
  double cdelt1, cdelt2;
  double crpix1, crpix2; /**< [pixel] */
  /** Coordinate value of reference pixel [m]. */
  double crval1, crval2; 
  
} PSF_Item;


/** Storage for the several PSFs available for a mirror system. */
typedef struct {
  /** Array of PSF_Items for different photon energies, off-axis
      angles, and azimuthal angles. */
  PSF_Item*** data;

  /** Number of different energies PSF images are available for. */
  int nenergies;
  /** Different energies PSF images are available for ([keV]). */
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


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Constructor for the PSF data structure. Reads PSF data from a FITS
    file with one or several image extensions. The file format is
    given by OGIP Calibration Memo CAL/GEN/92-027. The focal length of
    the telescope must be given. If the PSF image in the FITS file is
    defined with a pixel size in [arcsec], this value has to be
    converted to [m] using the knowledge about the focal length. */
PSF* newPSF(const char* const filename, const float focal_length,
	    int* const status);

/** Calculates the position on the detector, where a photon at given
    sky position with specified energy hits the detector according to
    the PSF data. The exact position is determined with a random
    number generator (randomization over one PSF pixel). The function
    return value is '1', if the photon hits the detector. If it does
    not fall onto the detector, the function returns '0'. The output
    detector position is stored in [m] in the first 2 parameters of
    the function. */
int get_psf_pos(/** Output: coordinates of the photon on the detector
		    ([m]). */
		struct Point2d* const position, 
		/** Incident photon. */
		const Photon photon, 
		/** Telescope information (pointing directions). */
		const struct Telescope telescope, 
		/** Focal length of the telescope [m]. Used for
		    off-axis angle calculations. */
		const float focal_length,
		/** Vignetting function. */
		const Vignetting* const vignetting, 
		/** PSF with data for different photon energies and
		    off-axis angles. */
		const PSF* const psf,
		int* const status);

/** Release the memory of the PSF storage. */
void destroyPSF(PSF** const psf);

/** Save the data contained in the PSF data struct to one or several
    image extensions in a FITS file following the OGIP standards for
    2-dimensional PSFs. */
int savePSFImage(const PSF* const psf, const char* const filename, 
		 int* const status);


#endif /* PSF_H */

