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



///////////////////////////////////////////////////////////////////////////////////


// This structure is used to store the PSF data for one particular 
// off-axis angle and one particular energy.
typedef struct {
  double **restrict data;   // pointer to PSF data array [x][y]

  double angle;             // off-axis angle of this particular PSF
  double energy;            // energy of this particular PSF

  // For each energy the PSFs are rescaled in such a way that the on-axis PSF 
  // is normalized to 1.  The scaling factor is stored here, and has to be used 
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


// Reads the PSF data file (containing FITS images) and
// stores this data in the PSF storage. 
PSF* get_psf(const char*, int* status);
//int get_psf_old(struct PSF_Store *store, const char filename[FILENAME_LENGTH]);


// Calculates the position on the detector, where a photon at given sky position 
// with specified energy hits the detector according to the PSF data and a random 
// number generator (randomization over one PSF pixel).
// Return value is '1', if the photon hits the detector. If it does not fall onto the
// detector, the function returns '0'.  The output detector position is stored 
// in [mu m] in the first 2 parameters of the function.
int get_psf_pos(struct Point2d*, struct Photon, struct Telescope, PSF*);


// Releases the memory for the PSF storage.
void free_psf(PSF*);


// Stores the PSF information in the PSF_Store to a fitsfile.
//int save_psf_to_fits(PSF *psf, const char filename[], int *status);


// creates the necessary parameters to create the table in the FITS file
static void psf_create_tbl_parameter(char *ftype[PSF_NFIELDS], 
				     char *fform[PSF_NFIELDS], 
				     char *funit[PSF_NFIELDS], int width);


// write PSF data to a FITS table
//int insert_psf_fitsrow(double angle, double energy, int x, int y, double *data, 
//		       long size, fitsfile *, long row);


// Save the data contained in the PSF storage to images in a FITS file
// according to the OGIP standards.
int save_psf_image(PSF*, const char *, int *status);


// This function reads one line of PSF data (i.e. a complete PSF image for one 
// particular off-axis angle and one particular photon energy) from a FITS table.
// The returned off-axis angle is given in [rad], the photon energy in [keV]. 
//int read_psf_fitsrow(double *angle, double *energy, int *x, int *y, double *data, 
//		     long size, fitsfile *fptr, long row);
//int read_psf_fitsrow(PSF_Item*, double *data, long size, fitsfile *fptr, 
//		     long row);


// Read the PSF images from a OGIP PSF FITS file.
//int read_psf_images(struct PSF_Store *, const char *, int *status);


#endif /* PSF_H */

