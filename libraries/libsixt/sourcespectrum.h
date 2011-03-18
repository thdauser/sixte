#ifndef SOURCESPECTRUM_H
#define SOURCESPECTRUM_H 1

#include "sixt.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Photon energy spectrum of an X-ray source. */
typedef struct {
  long nbins;
  float* emin;
  float* emax;

  /** Photon flux density [photons/s/cm^2/keV]. */
  float* flux;

  /** Accumulated photon rate distribution. The value of the last bin
      is the total photon rate distribution. In contrast to a common
      probability distribution, this array is not divided by the total
      photon rate, i.e., it is not normalized [integrated
      photons/s/bin] -> [photons/s]. */
  float* ratedistr;

  /** Filename the spectrum was loaded from. */
  char filename[MAXFILENAME];

} SourceSpectrum;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. */
SourceSpectrum* newSourceSpectrum(int* const status);

/** Destructor. */
void freeSourceSpectrum(SourceSpectrum** spec);

/** Return a random energy from a distribution according to the
    spectrum. */
float getRndSpectrumEnergy(const SourceSpectrum* const spec);


/** Load a spectrum from a file with the specified name. */
SourceSpectrum* loadXRaySpectrumFilename(const char* const filename,
					     int* const status);

/** Apply the instrument-specific ARF to a source spectrum in order to
    obtain a probability distribution from the given photon flux
    density. */
void applyARF2Spectrum(SourceSpectrum* const spec, 
		       const struct ARF* const arf, 
		       int* const status);

/** Determine the energy flux in the given energy band
    [erg/s/cm^2]. */
float getSpectralEnergyFlux(const SourceSpectrum* const spec,
			    const float emin, const float emax);

/** Determine the photon flux in the given energy band [photons/s]. */
float getSpectralPhotonRate(const SourceSpectrum* const spec,
			    const float emin, const float emax);


#endif /* SOURCESPECTRUM_H */
