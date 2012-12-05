#include "genpixgrid.h"


////////////////////////////////////////////////////////////////////
// Program Code
////////////////////////////////////////////////////////////////////


GenPixGrid* newGenPixGrid(int* const status) 
{
  // Allocate memory.
  GenPixGrid* grid=(GenPixGrid*)malloc(sizeof(GenPixGrid));
  if (NULL==grid) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for GenPixGrid failed");
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
  grid->rota  =0.;
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


void getGenDetAffectedPixel(const GenPixGrid* const grid, 
			    const double x,
			    const double y,
			    int* const xi,
			    int* const yi)
{
  // Calculate the integer indices.
  *xi=((int)((x+ (grid->xrpix-0.5)*grid->xdelt -grid->xrval)
	     /grid->xdelt +1.))-1;
  *yi=((int)((y+ (grid->yrpix-0.5)*grid->ydelt -grid->yrval)
	     /grid->ydelt +1.))-1;
  //                       |----|----> avoid (int)(-0.5) = 0


  // Check if this is a valid pixel.
  if ((*xi>=grid->xwidth) || (*xi<0) || (*yi>=grid->ywidth) || (*yi<0)) { 
    *xi=-1; 
    *yi=-1;
  } 

  // Check if the impact is located on the pixel border. 
  else if (grid->xborder>0.) {
    if ((grid->xrval+(*xi-grid->xrpix+1.5)*grid->xdelt-x<grid->xborder) ||
	(x-grid->xrval+(*xi-grid->xrpix+0.5)*grid->xdelt<grid->xborder)) {
      *xi=-1;
      *yi=-1;
    }
  }
  else if (grid->yborder>0.) {
    if ((grid->yrval+(*yi-grid->yrpix+1.5)*grid->ydelt-y<grid->yborder) ||
	(y-grid->yrval+(*yi-grid->yrpix+0.5)*grid->ydelt<grid->yborder)) {
      *xi=-1;
      *yi=-1;
    }
  }
}
