#ifndef POINTSOURCES_H
#define POINTSOURCES_H (1)

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fitsio.h"

#include "spectrum.h"
#include "vector.h"


// Maximum number of sources in the preselected source catalog:
#define MAX_N_POINTSOURCES 1000000


/** Contains all data to specify the properties of a point source. */
typedef struct {
  float ra, dec; /**< Right ascension and declination of the source [rad]. */
  float rate; /**< Average photon rate [photons/s]. */
  struct lightcurve_entry* lightcurve; /**< Pointer to source lightcurve. */
  struct Spectrum* spectrum; /**< Pointer to source spectrum. */
  //  struct PHA* pha_spectrum;   // source spectrum
  double t_last_photon; /** Time of last photon created for this source. */
} PointSource;

/** Contains all PointSources from one and the same FITS file. 
 * This data structure is generated from the sources sorted out of the 
 * PointSourceFile objects in the PointSourceFileCatalog. */
typedef struct {
  long nsources;
  PointSource* sources; /* nsources */
} PointSourceCatalog;



/** Gives information about a FITS file with point sources. */
typedef struct {
  fitsfile* fptr;

  long nrows; /**< Number of rows / PointSources in the FITS table. */

  /** Column names of the point source specific data columns in the FITS table. 
   * If a column is not available in the FITS file, the corresponding variable is 0.
   * If 'cspectrum' is 0, the default spectrum given by the header keyword 'SPECTRUM' 
   * should be used as source spectrum. */
  int cra, cdec, crate, cspectrum;
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

/** Constructor. */
PointSourceFile* get_PointSourceFile();
/** Enhanced constructor specifying a source file. */
PointSourceFile* get_PointSourceFile_fromFile(char* filename, int* status);
/** Destructor. */
void free_PointSourceFile(PointSourceFile*);

/** Selects point sources along the path of the telescope axis over the sky
 * and returns a PointSourceCatalog containing the individual PointSource objects. 
 * The selection is done based on a normal vector perpendicular to the scan plane and
 * the cos(sin?)-value of the maximum angle between the source direction and the scan plane.
 * For each selected source a Spectrum and a Lightcurve are assigned. */
PointSourceCatalog* get_PointSourceCatalog(PointSourceFileCatalog*, Vector normal_vector, 
					   const double max_align, struct Spectrum_Store,
					   int* status);
/** Destructor. */
void free_PointSourceCatalog(PointSourceCatalog* psc);


/** Reads a table from a FITS source catalogue. 
 * The right ascension and declination of the return PointSource structure are 
 * given in [rad]. */
int get_PointSourceTable_Row(PointSourceFile*, long row, PointSource*, int* status);

#endif /* POINTSOURCES_H */

