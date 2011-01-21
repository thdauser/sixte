#include "xraysource.h"


XRaySource* newXRaySource(int* const status)
{
  XRaySource* src = (XRaySource*)malloc(sizeof(XRaySource));
  CHECK_NULL(src, *status, 
	     "memory allocation for XRaySource failed");

  // Initalize pointers with NULL.
  src->t_next_photon = NULL;
  src->spectra       = NULL;

  // Initialize values.
  src->ra  = 0.;
  src->dec = 0.;
  src->pps = 0.;
  src->nspectra = 0;

  return(src);
}


void freeXRaySource(XRaySource** src)
{
  if (NULL!=*src) {
    if (NULL!=(*src)->t_next_photon) {
      free((*src)->t_next_photon);
    }
    if (NULL!=(*src)->spectra) {
      int ii;
      for (ii=0; ii<(*src)->nspectra; ii++) {
	freeXRaySourceSpectrum(&((*src)->spectra[ii]));
      }
    }
    free(*src);
    *src=NULL;
  }
}


LinkedListElement* getXRayPhotons(XRaySource* const src, 
				  const double t0, const double t1,
				  int* const status)
{
  // Time-ordered linked photon list.
  LinkedListElement* list=NULL;
  // Points to the last (NULL) pointer in the linked list.
  LinkedListElement** list_next = &list;

  // Photon arrival time.
  if (NULL==src->t_next_photon) {
    // There has been no photon for this particular source.
    src->t_next_photon=(double*)malloc(sizeof(double));
    CHECK_NULL(src->t_next_photon, *status,
	       "memory allocation for 't_next_photon' (double) failed");
    *(src->t_next_photon)= t0 + rndexp(1./src->pps);
  } else if (*(src->t_next_photon) < t0) {
    // The arrival time of the next photon was before the
    // requested time interval.
    *(src->t_next_photon)= t0 + rndexp(1./src->pps);
  }

  // Create new photons, as long as the requested time interval 
  // is not exceeded.
  while (*(src->t_next_photon) <= t1) {

    // Get a new photon object.
    Photon* ph = newPhoton(status);
    CHECK_STATUS_BREAK(*status);

    // Set the photon properties.
    ph->time = *(src->t_next_photon);
    ph->ra   = src->ra;
    ph->dec  = src->dec;

    // Determine the photon energy.
    ph->energy = getRndSpectrumEnergy(src->spectra[0]);

    // Append the photon at the end of the linked list.
    *list_next = newLinkedListElement(status);
    CHECK_STATUS_BREAK(*status);
    (*list_next)->el = ph;
    list_next = &((*list_next)->next);
  }
  CHECK_STATUS_RET(*status, list);

  // Return the linked list of photons.
  return(list);
}

