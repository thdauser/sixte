#include "genpixgrid.h"


GenPixGrid* newGenPixGrid(int* const status) 
{
  // Allocate memory.
  GenPixGrid* grid=(GenPixGrid*)malloc(sizeof(GenPixGrid));
  if (NULL==grid) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenPixGrid failed!\n", *status);
    return(grid);
  }

  // Initialize all pointers with NULL.

  return(grid);
}



void destroyGenPixGrid(GenPixGrid** const grid)
{
  if (NULL!=*grid) {
    free(*grid);
    *grid=NULL;
  }
}


