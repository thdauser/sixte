// The X-ray telescope simulation can handle several source spectra. 
// To each individual X-ray source a spectrum can be assigned. The 
// background spectrum is also handled in the same way as the source spectra.
// The simulation loads the spectra from PHA files containing the relative photon
// generation ratio per detector PHA channel. 
// That means in order to create an appropriate input 
// spectrum from a specified model, the energies have to be binned to the 
// EBOUNDS of the detector response.
// In principle the spectrum (i.e. the counts) can be normalized in an arbitrary way, as 
// the simulation normalizes the spectrum to 1 after reading it from the file, but it is
// strongly recommended to use only spectra that are normalized to 1.
// The absolute photon rate for a source has to be set via the 'PPS' field in the 
// source catalog.

#ifndef SPECTRUM_H
#define SPECTRUM_H 1

#include "sixt.h"
#include "fits_pha.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif


/** Data structure representing the specific spectrum of a source. */
typedef struct {
  /** Number of channels in the PHA spectrum. */
  long NumberChannels;

  /** Relative photon rates of the individual spectral bins. Array
      contains the photon rates per PHA channel. The rates respresent
      the relative photon rates, i.e., the total spectrum is
      normalized to 1. */
  float *rate; 
} Spectrum;


/** Storage for several different spectra. */
typedef struct {
  /** Array of the individual spectra in the storage. */
  Spectrum* spectrum; 
  /** Total number of spectra in the storage. */
  long nspectra; 
} SpectrumStore;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Load spectra from the files specified in a FITS header. The
    function scans the FITS header of the file given by the file
    pointer and searches for the NSPECTRA and SPECnnnn keywords. The
    specified source spectra are then loaded from the PHA files and
    store in the SpectrumStore. */
int loadSpectra(fitsfile*, SpectrumStore*);

/** Loads the specified PHA file containing a source spectrum. The
    data are read from the PHA file, the probability density
    (including normalization) of the relative photon rates in the
    individual channels is calculated and stored in the Spectrum data
    structure. */
int loadSpectrum(Spectrum*, char* filename);

/** Destructor for SpectrumStore. */
void freeSpectrumStore(SpectrumStore*);

/** Release previously allocated memory. */
void cleanupSpectrum(Spectrum*);


#endif  /* SPECTRUM_H */

