#ifndef ASTROSOURCES_TYPES_H
#define ASTROSOURCES_TYPES_H (1)

#include "vector.h"


// Define a structure which contains the necessary data for each source
struct source_cat_entry {
  struct vector r;                       // position of the source
  float rate;                            // photon rate
  struct lightcurve_entry* lightcurve;   // pointer to light curve
  struct Spectrum *spectrum;             // pointer to source spectrum

  double t_last_photon;      // time of last photon, which was created for this source 
};



#endif /* ASTROSOURCES_TYPES_H */

