#ifndef ASTROSOURCES_TYPES_H
#define ASTROSOURCES_TYPES_H (1)

#include "vector.h"
#include "fitsio.h"


// Maximum number of input source files:
#define MAX_N_POINTSOURCEFILES 5
// Maximum number of sources in the preselected source catalog:
#define MAX_N_POINTSOURCES 1000000



typedef struct {
  float ra, dec;  // right ascension and declination of the source [rad]
  float rate;     // photon rate
  struct lightcurve_entry* lightcurve; // pointer to source lightcurve
  struct Spectrum* spectrum;  // pointer to source spectrum REMOVE
  //  struct PHA* pha_spectrum;   // source spectrum
  double t_last_photon;       // time of last photon, which was created for this source
} PointSource;

typedef struct {
  PointSource* sources;
  long nsources;
} PointSourceCatalog;

typedef struct {
  fitsfile* files[MAX_N_POINTSOURCEFILES];
  int nfiles;
  int columns[MAX_N_POINTSOURCEFILES][3];
} PointSourceFiles;


#endif /* ASTROSOURCES_TYPES_H */

