// The X-ray telescope simulation can handle several source spectra. 
// To each individual X-ray source a spectrum can be assigned. The 
// background spectrum is also handled in the same way as the source spectra.
// The simulation loads the spectra from PHA files containing the relative photon
// generation ratio per detector PHA channel. 
// That means in order to create an appropriate input 
// spectrum from a specified model, the energies have to be binned to the 
// EBOUNDS of the detector response.
// The spectrum (i.e. the counts) can be normalized in an arbitrary way, as 
// the simulation normalizes the spectrum to 1 after reading it from the file. 
// The absolute photon rate for a source has to be set via the 'PPS' field in the 
// source catalog.

#ifndef SPECTRUM_H
#define SPECTRUM_H 1

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fitsio.h"
#include "headas.h"
#include "headas_error.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif

#include "fits_pha.h"
#include "sixt.h"


#define N_SPECTRA_FILES 1    // number of spectrum files (PHA files)


// Data structure representing the specific spectrum of a source.
struct Spectrum{
  long NumberChannels;

  float *rate; /**< Array of photon rates per PHA channel.
		* The rates respresent the relative photon rates.
                * The total spectrum is normalized to 1. */
};


// Storage for several different spectra.
struct Spectrum_Store{
  struct Spectrum* spectrum;
  long nspectra;

  //  struct PHA* pha_spectrum;
};



// Load the spectra from PHA files using the method "get_spectrum".
int get_spectra(struct Spectrum_Store *, long Nchannels, 
		char filenames[N_SPECTRA_FILES][FILENAME_LENGTH], int Nfiles);


// Loads the specified PHA file containing a source spectrum, 
// calculates the probability density (including normalization)
// and stores the spectrum in an array.
int get_spectrum(struct Spectrum *, long Nchannels, char filename[FILENAME_LENGTH]);

// Read a spectrum from a FITS PHA file (following the OGIP standards).
// The function uses routines from the HEAdas library 'libhdsp'.
// The spectrum is stored in the specified Spectrum_Store.
int assign_pha_spectrum(struct Spectrum_Store*, char* filename);

// Release memory of spectrum array
void free_spectra(struct Spectrum_Store *, long Nfiles);



#endif  /* SPECTRUM_H */

