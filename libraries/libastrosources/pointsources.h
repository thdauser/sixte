#ifndef POINTSOURCES_H
#define POINTSOURCES_H (1)

#include "sixt.h"
#include "spectrum.h"
#include "vector.h"
#include "linlightcurve.h"


/** Maximum number of sources in the preselected source catalog. */
#define MAX_N_POINTSOURCES 5000000


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** Contains all data to specify the properties of a point source. */
typedef struct {

  /** Right ascension and declination of the source [rad]. */
  float ra, dec; 

  /** Average photon rate [photons/s]. */
  float rate; 
  /** Type of the light curve particular for this PointSource. There
      are different possible values:

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
      SpectrumStore of the PointSourceCatalog. Warning: the
      spectrum_index starts at 1 (not at 0). */
  long spectrum_index; 
  /** Pointer to source spectrum. */
  Spectrum *spectrum; 

  /** Time of last photon created for this source. */
  double t_last_photon; 

} PointSource;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Sort the PointSourceList with the specified number of entries with
    respect to the requested coordinate axis using a quick sort
    algorithm. */
void quicksortPointSources(PointSource* list, long left, long right, int axis);


#endif /* POINTSOURCES_H */

