#include "pointsources.h"


long PointSourcesPartition(PointSource* list, 
			   long left, long right, 
			   long pivotIndex, int axis)
{
  Vector location;

  location = unit_vector(list[pivotIndex].ra, list[pivotIndex].dec);
  double pivotValue = getVectorDimensionValue(&location, axis);

  // Move pivot to end.
  PointSource buffer;
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



void quicksortPointSources(PointSource* list, long left, long right, int axis)
{
  if (right>left) {
    // select a pivot index //(e.g. pivotIndex := left+(right-left)/2)
    int pivotIndex = left+(right-left)/2;
    int pivotNewIndex = PointSourcesPartition(list, left, right, pivotIndex, axis);
    quicksortPointSources(list, left, pivotNewIndex-1, axis);
    quicksortPointSources(list, pivotNewIndex+1, right, axis);
  }
}



void freePointSource(PointSource ps)
{
  freeLinLightCurve(ps.lc);
}


