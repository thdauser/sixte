#ifndef POINTSOURCEFILE_H
#define POINTSOURCEFILE_H (1)

#include "sixt.h"
#include "vector.h"
#include "sourcelist.h"
#include "pointsources.h"


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** PointSourceFile containing information about the FITS file with
    the X-ray point sources. The FITS file contains a binary table
    with information about the different point sources: the position
    of the source, the reference photon rate of this source obtained
    from the spectral model and the instrument ARF, the number of the
    spectral model, and the number of the light curve model. The
    different spectral models are defined in the FITS header: the
    number of models is given by the keyword "NSPECTRA" and the
    filenames of the individual PHA files corresponding to the
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

  /** Number of rows / PointSources in the FITS table. */
  long nrows; 

  /** Column names of the point source specific data columns in the
      FITS table. If a column is not available in the FITS file, the
      corresponding variable is 0. */
  int cra, cdec, crate, cspectrum, clightcurve;

  /** SpectrumStore containing the spectra that are used in this
      PointSourceCatalog. */
  SpectrumStore spectrumstore;
} PointSourceFile;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor specifying a source file and the right HDU number to
    be loaded. */
PointSourceFile* get_PointSourceFile_fromFile(char* filename, int hdu, int* status);

/** Constructor for opening a PointSourceFile with the fitsfile
    pointer already pointing to the right HDU. */
PointSourceFile* get_PointSourceFile_fromHDU(fitsfile* fptr, int* status);

/** Destructor. */
void free_PointSourceFile(PointSourceFile*);

/** Reads a table from a FITS source catalogue. The line index of the
    requested line starts at 0 for the first line. The right ascension
    and declination of the returned PointSource structure are given in
    [rad]. */
int get_PointSourceTable_Row(PointSourceFile* psf, long row, 
			     PointSource* ps, int* status);

/** Read a SourceList from a given PointSourceFile. */
SourceList* readSourceListFromPointSourceFileHDU(fitsfile* file, 
						 long* nelements, 
						 int* status);


#endif /* POINTSOURCEFILE_H */

