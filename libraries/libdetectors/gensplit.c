#include "gensplit.h"


GenSplit* newGenSplit(int* const status) 
{
  // Allocate memory.
  GenSplit* split=(GenSplit*)malloc(sizeof(GenSplit));
  if (NULL==split) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenSplit failed!\n", *status);
    return(split);
  }

  // Initialize all pointers with NULL.

  // Set default values.
  split->type = GS_NONE;
  split->par1 = 0.;

  return(split);
}



void destroyGenSplit(GenSplit** const split)
{
  if (NULL!=*split) {
    free(*split);
    *split=NULL;
  }
}



void makeGenSplitEvents(const GenSplit* const split,
			const struct Point2d* const position,
			const float charge,
			const GenPixGrid* const grid,
			GenDetLine** const detline)
{
  // Which kind of split model has been selected?
  if (GS_NONE==split->type) {
    // No split events => all events are singles.

    // Determine the affected detector line and column.
    int line   = getGenDetAffectedLine  (grid, position->y);
    int column = getGenDetAffectedColumn(grid, position->x);

    // Check if the returned values are valid line and column indices.
    if ((0>column) || (0>line)) {
      return;
    }
    
    // TODO Introduce pile-up flag.
    
    // Add the charge (photon energy) to the affected pixel.
    addGenDetCharge2Pixel(detline[line], column, charge);
    
    // TODO Call the event trigger routine.
  
  } else {
    printf("Error: split model not supported!");
    exit(0);
  }
}

