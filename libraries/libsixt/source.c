#include "source.h"


Source* newSource(int* const status)
{
  Source* src = (Source*)malloc(sizeof(Source));
  CHECK_NULL(src, *status, 
	     "memory allocation for Source failed");

  // Initalize pointers with NULL.
  src->t_next_photon = NULL;
  src->ra            = 0.;
  src->dec           = 0.;
  src->extension     = 0.;
  src->row           = 0;

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
				     SimputCatalog* const simputcat,
				     const double t0, const double t1,
				     const double mjdref,
				     int* const status)
{
  // Time-ordered linked photon list.
  LinkedPhoListElement* list=NULL;
  // Points to the last (NULL) pointer in the linked list.
  LinkedPhoListElement** list_next = &list;

  // Load the source data from the SIMPUT catalog.
  SimputSource* simputsrc=loadCacheSimputSource(simputcat, src->row, status);
  CHECK_STATUS_RET(*status, list);

  // Photon arrival time.
  if (NULL==src->t_next_photon) {
    // There has been no photon for this particular source.
    src->t_next_photon=(double*)malloc(sizeof(double));
    CHECK_NULL(src->t_next_photon, *status,
	       "memory allocation for 't_next_photon' (double) failed");
    
    float rate=getSimputPhotonRate(simputsrc, t0, mjdref, status);
    CHECK_STATUS_RET(*status, list);
    if (0.==rate) {
      return(list);
    } else if (rate < 0.) {
      char msg[MAXMSG];
      sprintf(msg, "photon rate %e <= 0.", rate);
      SIXT_ERROR(msg);
      *status = EXIT_FAILURE;
      return(list);
    }
    *(src->t_next_photon) = t0 - 2./rate;

  } else if (*(src->t_next_photon) < t0) {
    float rate = getSimputPhotonRate(simputsrc, t0, mjdref, status);
    CHECK_STATUS_RET(*status, list);
    if (0.==rate) {
      return(list);
    } else if (rate <= 0.) {
      char msg[MAXMSG];
      sprintf(msg, "photon rate %e <= 0.", rate);
      SIXT_ERROR(msg);
      *status = EXIT_FAILURE;
      return(list);
    }
    *(src->t_next_photon) = t0 - 2./rate;
  }

  // Get a valid photon arrival time in the interval [t0,t1]
  while (*(src->t_next_photon) < t0) {
    // The arrival time of the next photon was before the
    // requested time interval.
    int failed=0;
    *(src->t_next_photon) = 
      getSimputPhotonTime(simputsrc, *(src->t_next_photon), mjdref, 
			  &failed, status);
    CHECK_STATUS_RET(*status, list);
    if (1==failed) return(list);
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
    getSimputPhotonCoord(simputsrc, &ph->ra, &ph->dec, status);
    CHECK_STATUS_RET(*status, list);
 
    // Copy the source identifiers.
    ph->src_id = simputsrc->src_id;

    // Set Photon ID to default value of 0. The proper value is
    // updated later, when the photon is inserted in the photon
    // list file.
    ph->ph_id = 0;

    // Determine the photon energy.
    ph->energy = 
      getSimputPhotonEnergy(simputsrc, *(src->t_next_photon), mjdref, status);
    CHECK_STATUS_BREAK(*status);

    // Determine the arrival time of the next (future) photon.
    int failed=0;
    *(src->t_next_photon) = 
      getSimputPhotonTime(simputsrc, *(src->t_next_photon), mjdref, 
			  &failed, status);
    CHECK_STATUS_BREAK(*status);
    if (1==failed) return(list);
  }
  CHECK_STATUS_RET(*status, list);

  // Return the linked list of photons.
  return(list);
}


static long SourcesPartition(Source* const list, 
			     const long left, const long right, 
			     const long pivotIndex, const int axis)
{
  Vector location = unit_vector(list[pivotIndex].ra, 
				list[pivotIndex].dec);
  double pivotValue = getVectorDimensionValue(&location, axis);

  // Move pivot to end.
  Source buffer;
  buffer = list[pivotIndex];
  list[pivotIndex] = list[right];
  list[right] = buffer;

  long storeIndex = left;
  long ii;
  for (ii=left; ii<right; ii++) { // left ≤ i < right  
    location = unit_vector(list[ii].ra, list[ii].dec);
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


