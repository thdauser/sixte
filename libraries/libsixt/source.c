#include "source.h"


Source* newSource(int* const status)
{
  Source* src = (Source*)malloc(sizeof(Source));
  CHECK_NULL(src, *status, 
	     "memory allocation for Source failed");

  // Initalize pointers with NULL.
  src->src           = NULL;
  src->t_next_photon = NULL;

  return(src);
}


void freeSource(Source** const src)
{
  if (NULL!=*src) {
    if (NULL!=(*src)->t_next_photon) {
      free((*src)->t_next_photon);
    }
    free(*src);
    *src=NULL;
  }
}


LinkedPhoListElement* getXRayPhotons(Source* const src, 
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
    *(src->t_next_photon) = getSimputPhotonTime(src->src, t0, status);
    CHECK_STATUS_RET(*status, list);

  } else if (*(src->t_next_photon) < t0) {
    // The arrival time of the next photon was before the
    // requested time interval.
    *(src->t_next_photon) = getSimputPhotonTime(src->src, t0, status);
    CHECK_STATUS_RET(*status, list);
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
    getSimputPhotonCoord(src->src, &ph->ra, &ph->dec, status);
    CHECK_STATUS_RET(*status, list);
 
    // Copy the source identifiers.
    ph->src_id = src->src->src_id;

    // Set Photon ID to default value of 0. The proper value is
    // updated later, when the photon is inserted in the photon
    // list file.
    ph->ph_id = 0;

    // Determine the photon energy.
    ph->energy = getSimputPhotonEnergy(src->src, status);
    CHECK_STATUS_BREAK(*status);

    // Determine the arrival time of the next (future) photon.
    *(src->t_next_photon) = 
      getSimputPhotonTime(src->src, *(src->t_next_photon), status);
    CHECK_STATUS_BREAK(*status);
  }
  CHECK_STATUS_RET(*status, list);

  // Return the linked list of photons.
  return(list);
}


static long SourcesPartition(Source* const list, 
			     const long left, const long right, 
			     const long pivotIndex, const int axis)
{
  Vector location = unit_vector(list[pivotIndex].src->ra, 
				list[pivotIndex].src->dec);
  double pivotValue = getVectorDimensionValue(&location, axis);

  // Move pivot to end.
  Source buffer;
  buffer = list[pivotIndex];
  list[pivotIndex] = list[right];
  list[right] = buffer;

  long storeIndex = left;
  long ii;
  for (ii=left; ii<right; ii++) { // left ≤ i < right  
    location = unit_vector(list[ii].src->ra, list[ii].src->dec);
    if (getVectorDimensionValue(&location, axis) <= pivotValue) {
      buffer = list[storeIndex];
      list[storeIndex] = list[ii];
      list[ii] = buffer;
      storeIndex++;
    }
  }

  // Move pivot to its final place
  buffer = list[storeIndex];
  list[storeIndex] = list[right];
  list[right] = buffer;

  return (storeIndex);
}


void quicksortSources(Source* const list, const long left, 
		      const long right, const int axis)
{
  if (right>left) {
    // select a pivot index //(e.g. pivotIndex := left+(right-left)/2)
    int pivotIndex = left+(right-left)/2;
    int pivotNewIndex = SourcesPartition(list, left, right, pivotIndex, axis);
    quicksortSources(list, left, pivotNewIndex-1, axis);
    quicksortSources(list, pivotNewIndex+1, right, axis);
  }
}


