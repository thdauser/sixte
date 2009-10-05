#ifndef LIGHTCURVE_H
#define LIGHTCURVE_H 1


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** Light curve giving the average source photon rate for a particular 
 * X-ray source. */
typedef struct {
  /** Number of time bins in the light curve. */
  long nbins;

  /** Width of each light curve bin ([s]). */
  double bin_width;

  /** Start time of the first light curve bin ([s]). */
  double t0;

  /** Light curve. Photon rate for each light curve bin ([photons/s]). */
  double* rate;

} LightCurve;


//////////////////////////////////////////////////////////////////////////
//   Function declarations
//////////////////////////////////////////////////////////////////////////


#endif /* LIGHTCURVE_H */
