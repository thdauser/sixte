#ifndef EXTENDEDSOURCES_H
#define EXTENDEDSOURCES_H (1)

#include "sixt.h"
#include "spectrum.h"
#include "vector.h"
#include "linlightcurve.h"


/** Maximum number of sources in the preselected source catalog. */
#define MAX_N_EXTENDEDSOURCES 1000000


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** Contains all data to specify the properties of an extended source. */
typedef struct {

  /** Right ascension and declination of the source [rad]. */
  float ra, dec; 

  /** Radial extension of the extended source [rad]. */
  float radius;

  /** Average photon rate [photons/s]. */
  float rate; 
  /** Type of the light curve particular for this
      ExtendedSource. There are different possible values:

      * T_LC_CONSTANT (=0) means a constant light curve.

      * T_LC_TIMMER_KOENIG (=-1) means a light curve with red noise
        according to Timmer & Koenig (1995).

      Positive values can be used to designate a particular FITS file
      containing a light curve. */
  long lc_type;
  /** Pointer to object with Piece-wise linear light curve for this
      X-ray source. */
  LinLightCurve* lc; 

  /** Index of the source spectrum within the SpectrumStore. This
      number represents the index of the source spectrum within the
      SpectrumStore of the ExtendedSourceCatalog. Warning: the
      spectrum_index starts at 1 (not at 0). */
  long spectrum_index; 
  /** Pointer to source spectrum. */
  Spectrum *spectrum; 

  /** Time of last photon created for this source. */
  double t_last_photon; 

} ExtendedSource;


/** Contains all ExtendedSources from one and the same FITS file. This
    data structure is generated from the sources sorted out of the
    ExtendedSourceFile. */
typedef struct {

  /** Number of ExtendedSource objects contained in the
      ExtendedSourceCatalog. */
  long nsources;

  /** Array containing the individual ExtendedSource objects. Length
      of the array is nsources. */
  ExtendedSource* sources;

} ExtendedSourceCatalog;


/** ExtendedSourceFile containing information about extended X-ray
    sources. Provides information about a FITS file with extended
    sources. The FITS file contains a binary table with information
    about the different extended sources: the position of the source,
    the characteristic source radius, the reference photon rate of
    this source obtained from the spectral model and the instrument
    ARF, the number of the spectral model, and the number of the light
    curve model. The different spectral models are defined in the FITS
    header: the number of models is given by the keyword "NSPECTRA"
    and the filenames of the individual PHA files corresponding to the
    spectral models by "SPECnnnn" with nnnn an integer number with 4
    digits starting from 0. The specification of the light curve
    models is similar: in the FITS header different light curve FITS
    files can be specified with the keywords "LCnnnnnn" with nnnnnn an
    integer number with 6 digits starting from 1. The column
    "lightcur" in the binary table can have integer values from -1 up
    to the maximum number of nnnnnn. The value -1 means that a
    red-noise light curve generated with the algorithm of Timmer &
    Koenig (1995) will be assigned to the source. The number 0 means
    that the source has a constant photon rate. For positive numbers
    the corresponding light curve model given in the FITS header will
    be assigned to the source. */
typedef struct {
  fitsfile* fptr;

  /** Number of rows / ExtendedSource's in the FITS table. */
  long nrows; 

  /** Column names of the extended source specific data columns in the
      FITS table. If a column is not available in the FITS file, the
      corresponding variable is 0. */
  int cra, cdec, cradius, crate, cspectrum, clightcurve;

  /** SpectrumStore containing the spectra that are used in this
      PointSourceCatalog. */
  SpectrumStore spectrumstore;
} PointSourceFile;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor returning a pointer to an empty ExtendedSourceFile object. */
ExtendedSourceFile* get_ExtendedSourceFile();

/** Enhanced constructor specifying a source file and the right HDU
    number to be loaded. */
ExtendedSourceFile* get_ExtendedSourceFile_fromFile(char* filename, int hdu, 
						    int* status);

/** Enhanced constructor for opening an ExtendedSourceFile with the
    fitsfile pointer already pointing to the right HDU. */
ExtendedSourceFile* get_ExtendedSourceFile_fromHDU(fitsfile* fptr, int* status);

/** Destructor. */
void free_ExtendedSourceFile(ExtendedSourceFile* esf);

/** Selects extended sources along the path of the telescope axis over
    the sky and returns an ExtendedSourceCatalog containing the
    individual ExtendedSource objects. The selection is done based on
    a normal vector perpendicular to the scan plane and the
    cos(sin?)-value of the maximum angle between the source direction
    and the scan plane. To each selected source a spectrum and a light
    curve are assigned. */
ExtendedSourceCatalog* getExtendedSourceCatalog(ExtendedSourceFile* esf, 
						Vector normal_vector, 
						const double max_align, 
						int* status);

/** Destructor. */
void free_ExtendedSourceCatalog(ExtendedSourceCatalog* esc);

/** Reads a table from a FITS extended source catalogue. The right
    ascension, the declination, and the radius of the returned
    ExtendedSource structure are given in [rad]. */
int get_ExtendedSourceTable_Row(ExtendedSourceFile* esf, long row, 
				ExtendedSource* es, int* status);


#endif /* EXTENDEDSOURCES_H */

