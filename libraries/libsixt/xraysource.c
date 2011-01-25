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


void freeXRaySource(XRaySource** const src)
{
  if (NULL!=*src) {
    if (NULL!=(*src)->t_next_photon) {
      free((*src)->t_next_photon);
    }
    if (NULL!=(*src)->spectra) {
      free((*src)->spectra);
      // NOTE: the spectra must not be free'd, because this
      // is already done in the destructor of the XRaySourceCatalog.
    }
    free(*src);
    *src=NULL;
  }
}


LinkedPhoListElement* getXRayPhotons(XRaySource* const src, 
				     const double t0, const double t1,
				     int* const status)
{
  // Time-ordered linked photon list.
  LinkedPhoListElement* list=NULL;
  // Points to the last (NULL) pointer in the linked list.
  LinkedPhoListElement** list_next = &list;

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

    // Append a new entry at the end of the linked list.
    *list_next = newLinkedPhoListElement(status);
    CHECK_STATUS_BREAK(*status);
    Photon* ph = &((*list_next)->photon);
    list_next =  &((*list_next)->next);

    // Set the photon properties.
    ph->time = *(src->t_next_photon);
    ph->ra   = src->ra;
    ph->dec  = src->dec;

    // Determine the photon energy.
    ph->energy = getRndSpectrumEnergy(src->spectra[0]);
    
    // Determine the arrival time of the next (future) photon.
    *(src->t_next_photon) += rndexp(1./src->pps);
  }
  CHECK_STATUS_RET(*status, list);

  // Return the linked list of photons.
  return(list);
}


static long XRaySourcesPartition(XRaySource* const list, 
				 const long left, const long right, 
				 const long pivotIndex, const int axis)
{
  Vector location = unit_vector(list[pivotIndex].ra, list[pivotIndex].dec);
  double pivotValue = getVectorDimensionValue(&location, axis);

  // Move pivot to end.
  XRaySource buffer;
  buffer = list[pivotIndex];
  list[pivotIndex] = list[right];
  list[right] = buffer;

  long storeIndex = left;
  long i;
  for (i=left; i<right; i++) { // left ≤ i < right  
    location = unit_vector(list[i].ra, list[i].dec);
    if (getVectorDimensionValue(&location, axis) <= pivotValue) {
      buffer = list[storeIndex];
      list[storeIndex] = list[i];
      list[i] = buffer;
      storeIndex++;
    }
  }

  // Move pivot to its final place
  buffer = list[storeIndex];
  list[storeIndex] = list[right];
  list[right] = buffer;

  return (storeIndex);
}


void quicksortXRaySources(XRaySource* const list, const long left, 
			  const long right, const int axis)
{
  if (right>left) {
    // select a pivot index //(e.g. pivotIndex := left+(right-left)/2)
    int pivotIndex = left+(right-left)/2;
    int pivotNewIndex = XRaySourcesPartition(list, left, right, pivotIndex, axis);
    quicksortXRaySources(list, left, pivotNewIndex-1, axis);
    quicksortXRaySources(list, pivotNewIndex+1, right, axis);
  }
}


