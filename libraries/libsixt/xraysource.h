#ifndef XRAYSOURCE_H
#define XRAYSOURCE_H 1

#include "sixt.h"

#include "linkedpholist.h"
#include "photon.h"
#include "xraysourcespectrum.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Photon energy spectrum of an X-ray source. */
typedef struct {

  float ra, dec;

  /** Photon rate [photons/s]. */
  float pps;

  /** Time of the emission of the last photon. */
  double* t_next_photon;

  /** Source spectrum / spectra. */
  int nspectra;
  XRaySourceSpectrum** spectra;

  /* TODO
  int nimages;
  XRaySourceImage** images;

  XRayLightCurve* lc;
  */

} XRaySource;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. */
XRaySource* newXRaySource(int* const status);

/** Destructor. */
void freeXRaySource(XRaySource** const src);

/** Create photons for a particular source in the specified time
    interval. */
LinkedPhoListElement* getXRayPhotons(XRaySource* const src, 
				     const double t0, const double t1,
				     int* const status);

/** Sort the list of XRaySource objects with the specified number of
    entries with respect to the requested coordinate axis using a
    quick sort algorithm. */
void quicksortXRaySources(XRaySource* const list, const long left, 
			  const long right, const int axis);


#endif /* XRAYSOURCE_H */
