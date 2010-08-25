#include "gendet.h"


GenDet* newGenDet(int* status) {
  GenDet* det=NULL;

  // Allocate memory.
  det=(GenDet*)malloc(sizeof(GenDet));
  if (NULL==det) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenDet failed!\n", *status);
    return(det);
  }

  // Initialize all pointers with NULL.
  det->line=NULL;
  det->rmf =NULL;

  // TODO Read in the XML definition of the detector.

  // Allocate memory for the pixels.
  det->line=(GenDetLine**)malloc(384*sizeof(GenDetLine*));
  if (NULL==det->line) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenDet failed!\n", *status);
    return(det);
  }
  int i;
  for (i=0; i<384; i++) {
    det->line[i] = newGenDetLine(384, status);
    if (EXIT_SUCCESS!=*status) return(det);
  }

  return(det);
}



void destroyGenDet(GenDet** gendet);
