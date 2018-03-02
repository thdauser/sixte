/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
*/

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
  grid->xwidth=INT_MAX;
  grid->ywidth=INT_MAX;
  grid->xrpix =NAN;
  grid->yrpix =NAN;
  grid->xrval =NAN;
  grid->yrval =NAN;
  grid->xdelt =NAN;
  grid->ydelt =NAN;
  grid->rota  =NAN;
  grid->xborder=0.; // dead area: 0 by default
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
			    int* const yi,
			    double* const xp,
			    double* const yp)
{
  // Calculate the distance to the reference point [m].
  double xd=x-grid->xrval;
  double yd=y-grid->yrval;

  // Rotate around the reference point [m].
  double cosrota=cos(grid->rota);
  double sinrota=sin(grid->rota);
  double xr= xd*cosrota +yd*sinrota;
  double yr=-xd*sinrota +yd*cosrota;

  // Calculate the real valued pixel indices.
  double xb=xr/grid->xdelt + (grid->xrpix+0.5);
  double yb=yr/grid->ydelt + (grid->yrpix+0.5);

  // Calculate the integer pixel indices.
  *xi=((int)(xb +1.))-1;
  *yi=((int)(yb +1.))-1;
  //             |----|---->  avoid (int)(-0.5) = 0


  // Check if this is a valid pixel.
  if ((*xi>=grid->xwidth) || (*xi<0) || (*yi>=grid->ywidth) || (*yi<0)) { 
    *xi=-1; 
    *yi=-1;
    return;
  } else { 
    // Check if the impact is located on one of the pixel borders
    if (grid->xborder>0. || grid->yborder>0.) {
      if ((grid->xrval+(*xi-grid->xrpix+0.5)*grid->xdelt-x<grid->xborder) ||
	  (x-grid->xrval+(*xi-grid->xrpix-0.5)*grid->xdelt<grid->xborder)) {
	*xi=-1;
	*yi=-1;
	return;
      }
      if ((grid->yrval+(*yi-grid->yrpix+0.5)*grid->ydelt-y<grid->yborder) ||
	  (y-grid->yrval+(*yi-grid->yrpix-0.5)*grid->ydelt<grid->yborder)) {
	*xi=-1;
	*yi=-1;
	return;
      }
    }
  }

  // Calculate the relative position with respect to the left 
  // and lower pixel boundaries.
  *xp=xb-1.0*(*xi);
  *yp=yb-1.0*(*yi);
}

