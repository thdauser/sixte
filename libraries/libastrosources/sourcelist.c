#include "sourcelist.h"


long SourceListPartition(SourceList* list, 
			 long left, long right, 
			 long pivotIndex, int axis)
{
  double pivotValue = 
    getVectorDimensionValue(&list[pivotIndex].location, axis);

  // Move pivot to end.
  struct SourceListEntry buffer;
  buffer = list[pivotIndex];
  list[pivotIndex] = list[right];
  list[right] = buffer;

  long storeIndex = left;
  long i;
  for (i=left; i<right; i++) { // left ≤ i < right  
    if (getVectorDimensionValue(&list[i].location, axis) <= pivotValue) {
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


void quicksortSourceList(SourceList* list, long left, long right, int axis)
{
  if (right>left) {
    // select a pivot index //(e.g. pivotIndex := left+(right-left)/2)
    int pivotIndex = left+(right-left)/2;
    int pivotNewIndex = SourceListPartition(list, left, right, pivotIndex, axis);
    quicksortSourceList(list, left, pivotNewIndex-1, axis);
    quicksortSourceList(list, pivotNewIndex+1, right, axis);
  }
}


