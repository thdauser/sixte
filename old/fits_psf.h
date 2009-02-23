#ifndef FITS_PSF_H
#define FITS_PSF_H 1


#include "fitsio.h"

//#include "psf.h"

#define PSF_N_ANGLES 7      // number of different off-axis angles
#define PSF_N_ENERGIES 3    // number of different energies
#define PSF_NFIELDS 5       // number of fields in the table (= number of columns)

#define MAXMSG 256

/*
// Store the information about the index - energy and index - off-axis-angle relations in the PSF data structure,
// in order to be able look up, which energy/angle corresponds to which index:
const double angles[PSF_N_ANGLES] = {0.,                  //  [floating point pixel]
				     5./60.,
				     10./60.,
				     15./60.,
				     20./60.,
				     25./60.,
				     30./60. };

const double energies[PSF_N_ENERGIES] = {1., 4., 7.};     // energies in [keV]
*/



// creates the necessary parameters to create the table in the FITS file
void psf_create_tbl_parameter(
			      char *ftype[PSF_NFIELDS], char *fform[PSF_NFIELDS], char *funit[PSF_NFIELDS],
			      int width);

// write PSF data to a FITS table
int insert_psf_fitsrow(double angle, double energy, int x, int y, double *data, long size, fitsfile *output_fptr, long row);

// read PSF data from a FITS table
int read_psf_fitsrow(double *angle, double *energy, int *x, int *y, double *data, long size, fitsfile *fptr, long row);


#endif

