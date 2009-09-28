#ifndef POINTSOURCES_H
#define POINTSOURCES_H (1)

#include "sixt.h"
#include "spectrum.h"
#include "vector.h"


// Maximum number of sources in the preselected source catalog:
#define MAX_N_POINTSOURCES 1000000


/** Contains all data to specify the properties of a point source. */
typedef struct {
  float ra, dec; /**< Right ascension and declination of the source [rad]. */
  float rate; /**< Average photon rate [photons/s]. */
  struct lightcurve_entry* lightcurve; /**< Pointer to source lightcurve. */

  Spectrum *spectrum; /**< Pointer to source spectrum. */

  /** Index of the source spectrum within the SpectrumStore.
   * This number represents the index of the source spectrum within the SpectrumStore
   * of the PointSourceCatalog. 
   * Warning: the spectrum_index starts  at 1 (not at 0). */
  long spectrum_index; 
			  

  double t_last_photon; /** Time of last photon created for this source. */
} PointSource;

/** Contains all PointSources from one and the same FITS file. 
 * This data structure is generated from the sources sorted out of the 
 * PointSourceFile objects in the PointSourceFileCatalog. */
typedef struct {
  /** Number of PointSource objects contained in the PointSourceCatalog. */
  long nsources;
  /** Array containing the individual PointSource objects. 
   * Length of the array is nsources. */
  PointSource* sources;
} PointSourceCatalog;


/** Gives information about a FITS file with point sources. */
typedef struct {
  fitsfile* fptr;

  long nrows; /**< Number of rows / PointSources in the FITS table. */

  /** Column names of the point source specific data columns in the FITS table. 
   * If a column is not available in the FITS file, the corresponding variable is 0. */
  int cra, cdec, crate, cspectrum;

  /** SpectrumStore containing the spectra that are used in this PointSourceCatalog. */
  SpectrumStore spectrumstore;
} PointSourceFile;

/** Collection of point source FITS files.
 * Used as input for the point source selection algorithm. */
typedef struct {
  int nfiles;
  PointSourceFile** files;
} PointSourceFileCatalog;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. */
PointSourceFileCatalog* get_PointSourceFileCatalog();
/** Destructor. */
void free_PointSourceFileCatalog(PointSourceFileCatalog*);

/** Constructor returning a pointer to an empty PointSourceFile object. */
PointSourceFile* get_PointSourceFile();
/** Enhanced constructor specifying a source file to be loaded. */
PointSourceFile* get_PointSourceFile_fromFile(char* filename, int* status);
/** Destructor. */
void free_PointSourceFile(PointSourceFile*);

/** Selects point sources along the path of the telescope axis over the sky
 * and returns a PointSourceCatalog containing the individual PointSource objects. 
 * The selection is done based on a normal vector perpendicular to the scan plane and
 * the cos(sin?)-value of the maximum angle between the source direction and the scan plane.
 * For each selected source a Spectrum and a Lightcurve are assigned. */
PointSourceCatalog* get_PointSourceCatalog(PointSourceFileCatalog*, Vector normal_vector, 
					   const double max_align, int* status);
/** Destructor. */
void free_PointSourceCatalog(PointSourceCatalog* psc);

/** Reads a table from a FITS source catalogue. 
 * The right ascension and declination of the return PointSource structure are 
 * given in [rad]. */
int get_PointSourceTable_Row(PointSourceFile*, long row, PointSource*, int* status);

#endif /* POINTSOURCES_H */

