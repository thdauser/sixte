#include "genpixgrid.h"


////////////////////////////////////////////////////////////////////
// Static function declarations
////////////////////////////////////////////////////////////////////


/** Return the index of the bin affected by the specified
    x-position. The bin grid is defined by the following WCS compliant
    values: reference pixel (rpix), reference value (rval), and pixel
    delta (delt). If the specified x-value is outside the bins, the
    function return value is -1. */
static inline int getAffectedIndex(const double x, const float rpix, 
				   const float rval, const float delt, 
				   const int width);


////////////////////////////////////////////////////////////////////
// Program Code
////////////////////////////////////////////////////////////////////


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



static inline int getAffectedIndex(const double x, const float rpix, 
				   const float rval, const float delt, 
				   const int width)
{
  int index = ((int)((x+ (rpix-0.5)*delt -rval)/delt +1.))-1;
  //                  avoid (int)(-0.5) = 0     <-----|----|
  if (index>=width) { index=-1; }
  return(index);
}



inline int getGenDetAffectedLine(const GenPixGrid* const grid, const double y)
{
  return(getAffectedIndex(y, grid->yrpix, grid->yrval, grid->ydelt, grid->ywidth));
}



inline int getGenDetAffectedColumn(const GenPixGrid* const grid, const double x)
{
  return(getAffectedIndex(x, grid->xrpix, grid->xrval, grid->xdelt, grid->xwidth));
}


