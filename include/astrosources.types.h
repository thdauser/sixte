#ifndef ASTROSOURCES_TYPES_H
#define ASTROSOURCES_TYPES_H (1)

#include "vector.h"
#include "fitsio.h"


// Maximum number of input source files:
#define MAX_N_POINTSOURCEFILES 5
#define MAX_NSOURCEFILES 5
// Maximum number of sources in the preselected source catalog:
#define MAX_NSOURCES_PRE 1000000




// Define a structure which contains the necessary data for each source
struct source_cat_entry {
  double ra, dec;   // right ascension and declination of the source position
  struct vector r;  // REMOVE            // position of the source
  
  float rate;                            // photon rate
  struct lightcurve_entry* lightcurve;   // pointer to light curve
  struct Spectrum *spectrum;             // pointer to source spectrum

  double t_last_photon;      // time of last photon, which was created for this source 
};



typedef struct {
  float ra, dec;  // right ascension and declination of the source
  float rate;     // photon rate
  struct lightcurve_entry* lightcurve; // pointer to source lightcurve
  struct Spectrum* spectrum;  // pointer to source spectrum
  double t_last_photon;       // time of last photon, which was created for this source
} PointSource;

typedef struct {
  PointSource* sources;
} PointSourceCatalog;

typedef struct {
  fitsfile* files[MAX_N_POINTSOURCEFILES];
  int columns[MAX_N_POINTSOURCEFILES][3];
} PointSourceFiles;


#endif /* ASTROSOURCES_TYPES_H */

