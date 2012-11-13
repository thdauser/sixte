#include "genpixgrid.h"


////////////////////////////////////////////////////////////////////
// Static function declarations
////////////////////////////////////////////////////////////////////


/** Return the index of the bin affected by the specified
    x-position. The bin grid is defined by the following WCS compliant
    values: reference pixel (rpix), reference value (rval), and pixel
    delta (delt). If the specified x-value is outside the bins, the
    function return value is -1. If the photon impact lies on the
    pixel border the return valus is -1. */
static inline int getAffectedIndex(const double x, const float rpix, 
				   const float rval, const float delt, 
				   const float border, const int width);


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

  // Initialize values.
  grid->xwidth=0;
  grid->ywidth=0;
  grid->xrpix =0.;
  grid->yrpix =0.;
  grid->xrval =0.;
  grid->yrval =0.;
  grid->xdelt =0.;
  grid->ydelt =0.;
  grid->xborder=0.;
  grid->yborder=0.;

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
				   const float border, const int width)
{
  int index = ((int)((x+ (rpix-0.5)*delt -rval)/delt +1.))-1;
  //                  avoid (int)(-0.5) = 0     <-----|----|
  if (index>=width) { 
    index=-1; 
  } else if (border>0.) {
    // Check if the impact is located on the pixel border. 
    if ((rval+(index-rpix+1.5)*delt-x<border) ||
	(x-rval+(index-rpix+0.5)*delt<border)) {
      index=-1;
    }
  }
  return(index);
}


int getGenDetAffectedLine(const GenPixGrid* const grid, const double y)
{
  return(getAffectedIndex(y, grid->yrpix, grid->yrval, grid->ydelt, 
			  grid->yborder, grid->ywidth));
}


int getGenDetAffectedColumn(const GenPixGrid* const grid, const double x)
{
  return(getAffectedIndex(x, grid->xrpix, grid->xrval, grid->xdelt, 
			  grid->xborder, grid->xwidth));
}

