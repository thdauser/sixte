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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "squarepixels.h"


SquarePixels* newSquarePixels(struct SquarePixelsParameters* spp, int* const status)
{
  //Memory-allocation
  SquarePixels* sp=(SquarePixels*)malloc(sizeof(SquarePixels));
  if (NULL==sp) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for SquarePixels "
		   "data structure!\n", *status);
    return(sp);
  }
  //Set initial values.
  sp->array=NULL;
  sp->line2readout=NULL;

  //Set the array dimensions via SquarePixelsParameters:
  sp->xwidth = spp->xwidth;
  sp->ywidth = spp->ywidth;
  sp->xoffset = spp->xwidth/2;
  sp->yoffset = spp->ywidth/2;
  sp->xpixelwidth = spp->xpixelwidth;
  sp->ypixelwidth = spp->ypixelwidth;
  sp->DCU_length = spp->DCU_length;
  if(spp->DCU_length!=0.){
  sp->DCU_gap = spp->DCU_gap;
  sp->DCA_gap = spp->DCA_gap;
  }

  //Get the memory for the pixels(xwidth x ywidth):
  //Rows:size of one pixel * number of pixels in x-direction
  sp->array = (SquarePixel**)malloc(sp->xwidth*sizeof(SquarePixel*));
  if (NULL==sp->array) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: memory allocation for pixel array failed!\n",
		   *status);
    return(sp);
  } else {
    int count;
    //Columns:size of one pixel * number of pixels in y-direction,
    //for each column
    for(count=0; count<sp->xwidth; count++) {
      sp->array[count] = (SquarePixel*)malloc(sp->ywidth*sizeof(SquarePixel));
      if (NULL==sp->array[count]) {
	*status = EXIT_FAILURE;
	HD_ERROR_THROW("Error: memory allocation for pixel array failed!\n",
		       *status);
	return(sp);
      }
    }
  }

  // Get the memory for the read-out flags of the individual detector lines.
  sp->line2readout = (int*)malloc(sp->ywidth*sizeof(int));
  if (NULL==sp->line2readout) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: memory allocation for pixel array failed!\n",
		   *status);
    return(sp);
  }

  //Clear the pixels.
  //First row:clears all columns, next row: all columns, ...
  clearSquarePixels(sp);

  return(sp);
}


void clearLineSquarePixels(SquarePixels* const sp, const int line)
{
  int x;
  for (x=0; x<sp->xwidth; x++) {
    sp->array[x][line].charge = 0.;
    sp->array[x][line].valid_flag = 0;
  }
  sp->line2readout[line] = 0;
}


void clearSquarePixels(SquarePixels* const sp)
{
  int line;
  for (line=0; line<sp->ywidth; line++) {
    clearLineSquarePixels(sp, line);
  }
}


void destroySquarePixels(SquarePixels** const sp)
{
  if (NULL!=*sp) {
    if (NULL!=(*sp)->array) {
      int count;
      for (count=0; count<(*sp)->xwidth; count++) {
	if (NULL!=(*sp)->array[count]) {
	  free((*sp)->array[count]);
	}
      }
      free((*sp)->array);
      (*sp)->array=NULL;
    }

    if (NULL!=(*sp)->line2readout) {
      free((*sp)->line2readout);
      (*sp)->line2readout=NULL;
    }

    free(*sp);
  }
}


/** Determines the minimum distance value out of an array with 4
    entries and returns the corresponding index. */
static inline int getMinimumDistance(double array[])
{
  int count, index=0;
  double minimum=array[0];

  for (count=1; count<4; count++) {
    if ( (minimum < 0.) ||
	 ((array[count]<=minimum)&&(array[count]>=0.)) ) {
      minimum = array[count];
      index = count;
    }
  }

  return(index);
}


int getSquarePixelsGaussianSplits(SquarePixels* sp, GaussianChargeCloud* gcc,
				  struct Point2d position,
				  int* x, int* y, double* fraction)
{
  int npixels = 0; // Number of affected pixels.

  // The following array entries are used to transform between
  // different array indices.
  int xe[4] = {1, 0,-1, 0};
  int ye[4] = {0, 1, 0,-1};

  // Calculate pixel indices (integer) of central affected pixel:
  x[0] = (int)(position.x/sp->xpixelwidth + 0.5*sp->xwidth +1.)-1;
  y[0] = (int)(position.y/sp->ypixelwidth + 0.5*sp->ywidth +1.)-1;

  // Check if the main event lies inside the detector.
  if ((x[0]<0) || (x[0]>=sp->xwidth) || (y[0]<0) || (y[0]>=sp->ywidth)) {
    x[0] = INVALID_PIXEL;
    y[0] = INVALID_PIXEL;
    return(0);
  }

  // If charge cloud size is 0, i.e. no splits are created.
  if (gcc->ccsize < 1.e-20) {
    // Only single events are possible!
    npixels = 1;
    fraction[0] = 1.;

  } else {
    // Calculate the distances from the event to the borders of the
    // surrounding pixel (in [m]):
    double distances[4] = {
      // distance to right pixel edge
      (x[0]-sp->xoffset+1)*sp->xpixelwidth - position.x,
      // distance to upper edge
      (y[0]-sp->yoffset+1)*sp->ypixelwidth - position.y,
      // distance to left pixel edge
      position.x - (x[0]-sp->xoffset)*sp->xpixelwidth,
      // distance to lower edge
      position.y - (y[0]-sp->yoffset)*sp->ypixelwidth
    };

    int mindist = getMinimumDistance(distances);
    if (distances[mindist] < gcc->ccsize) {
      // Not a single event!
      x[1] = x[0] + xe[mindist];
      y[1] = y[0] + ye[mindist];

      double mindistgauss = gaussint(distances[mindist]/gcc->ccsigma);

      double minimum = distances[mindist];
      // search for the next to minimum distance to an edge
      distances[mindist] = -1.;
      int secmindist = getMinimumDistance(distances);
      distances[mindist] = minimum;

      if (distances[secmindist] < gcc->ccsize) {
	// Quadruple!
	npixels = 4;

	x[2] = x[0] + xe[secmindist];
	y[2] = y[0] + ye[secmindist];
	x[3] = x[1] + xe[secmindist];
	y[3] = y[1] + ye[secmindist];

	// Calculate the different charge fractions in the 4 affected pixels.
	double secmindistgauss = gaussint(distances[secmindist]/gcc->ccsigma);
	fraction[0] = (1.-mindistgauss)*(1.-secmindistgauss);
	fraction[1] =     mindistgauss *(1.-secmindistgauss);
	fraction[2] = (1.-mindistgauss)*    secmindistgauss ;
	fraction[3] =     mindistgauss *    secmindistgauss ;

      } else {
	// Double!
	npixels = 2;

	fraction[0] = 1. - mindistgauss;
	fraction[1] =      mindistgauss;

      } // END of double or Quadruple

    } else {
      // Single event!
      npixels = 1;
      fraction[0] = 1.;

    } // END of check for single event
  } // END of check for charge cloud size equals 0.

  // Check whether all split partners lie inside the detector:
  int count;
  for(count=1; count<npixels; count++) {
    if ((x[count]<0) || (x[count]>=sp->xwidth) ||
	(y[count]<0) || (y[count]>=sp->ywidth)) {
      x[count] = INVALID_PIXEL;
      y[count] = INVALID_PIXEL;
    }
  }

  // Return the number of affected pixels, or 0, if the position lies outside
  // the detector array.
  return(npixels);
}



int getSquarePixelsExponentialSplits(SquarePixels* sp, ExponentialChargeCloud* ecc,
				     struct Point2d position,
				     int* x, int* y, double* fraction)
{
  // None-Gaussian, exponential charge cloud model
  // (concept proposed by Konrad Dennerl).
  int npixels = 4; // Number of affected pixels.

  // The following array entries are used to transform between
  // different array indices.
  int xe[4] = {1, 0,-1, 0};
  int ye[4] = {0, 1, 0,-1};

  // Determine the affected pixels:

  // Calculate pixel indices (integer) of central affected pixel:
  x[0] = (int)(position.x/sp->xpixelwidth + (double)(sp->xoffset)+1.) -1;
  y[0] = (int)(position.y/sp->ypixelwidth + (double)(sp->yoffset)+1.) -1;

  // Check if the main event lies inside the detector.
  if ((x[0]<0) || (x[0]>=sp->xwidth) || (y[0]<0) || (y[0]>=sp->ywidth)) {
    x[0] = INVALID_PIXEL;
    y[0] = INVALID_PIXEL;
    return(0);
  }

  // Calculate the distances from the exact event position to the borders of the
  // surrounding pixel (in units [fraction of a pixel border]):
  double distances[4] = {
    // distance to right pixel edge
    (x[0]-sp->xoffset+1)*1. - position.x/sp->xpixelwidth,
    // distance to upper edge
    (y[0]-sp->yoffset+1)*1. - position.y/sp->ypixelwidth,
    // distance to left pixel edge
    position.x/sp->xpixelwidth - (x[0]-sp->xoffset)*1.,
    // distance to lower edge
    position.y/sp->ypixelwidth - (y[0]-sp->yoffset)*1.
  };

  // Search for the minimum distance to the edges.
  int mindist = getMinimumDistance(distances);
  x[1] = x[0] + xe[mindist];
  y[1] = y[0] + ye[mindist];

  // Search for the next to minimum distance to the edges.
  double minimum = distances[mindist];
  distances[mindist] = -1.;
  int secmindist = getMinimumDistance(distances);
  distances[mindist] = minimum;
  // Pixel coordinates of the 3rd and 4th split partner.
  x[2] = x[0] + xe[secmindist];
  y[2] = y[0] + ye[secmindist];
  x[3] = x[1] + xe[secmindist];
  y[3] = y[1] + ye[secmindist];

  // Now we know the affected pixels and can determine the
  // charge fractions according to the model exp(-(r/0.355)^2).
  // Remember that the array distances[] contains the distances
  // to the pixel borders, whereas here we need the distances from
  // the pixel center for the parameter r.
  // The value 0.355 is given by the parameter ecc->parameter.
  fraction[0] = exp(-(pow(0.5-distances[mindist],2.)+pow(0.5-distances[secmindist],2.))/
		    pow(ecc->parameter,2.));
  fraction[1] = exp(-(pow(0.5+distances[mindist],2.)+pow(0.5-distances[secmindist],2.))/
		    pow(ecc->parameter,2.));
  fraction[2] = exp(-(pow(0.5-distances[mindist],2.)+pow(0.5+distances[secmindist],2.))/
		    pow(ecc->parameter,2.));
  fraction[3] = exp(-(pow(0.5+distances[mindist],2.)+pow(0.5+distances[secmindist],2.))/
		    pow(ecc->parameter,2.));
  // Normalization to 1.
  double sum = fraction[0]+fraction[1]+fraction[2]+fraction[3];
  fraction[0] /= sum;
  fraction[1] /= sum;
  fraction[2] /= sum;
  fraction[3] /= sum;

  // Check whether all split partners lie inside the detector:
  int count;
  for (count=1; count<npixels; count++) {
    if ((x[count]<0) || (x[count]>=sp->xwidth) ||
	(y[count]<0) || (y[count]>=sp->ywidth)) {
      x[count] = INVALID_PIXEL;
      y[count] = INVALID_PIXEL;
    }
  }

  // Return the number of affected pixels.
  return(npixels);
}


int getSquarePixel(SquarePixels* sp, struct Point2d position, int* x, int* y)
{
  if(sp->DCU_length!=0.){//detector with gaps

  //position.x and .y is with reference to the back left corner of the detector [m]
  //(if pos. x-dir is from back to front and pos- y-dir from left to right)
  //first: check, whether impact position lies within one of the DCU's
  //second: determine position in terms of pixels
  //third: ensure that there will be no error with the type-cast
  //(int rounds to lower value)

  double DCA_length=(2.*sp->DCU_length)+sp->DCU_gap;
  double DCA=DCA_length+sp->DCA_gap;

  int ratio_x=(int)((position.x)/DCA);
  int ratio_y=(int)((position.y)/DCA);

  if(position.x < (ratio_x*DCA+sp->DCU_length))
    {*x = (int)(position.x/sp->xpixelwidth +1.)-1;}
     else{if(position.x > (ratio_x*DCA+sp->DCU_length+sp->DCU_gap)
			   && position.x < (ratio_x*DCA+DCA_length))
	 {*x = (int)(position.x/sp->xpixelwidth +1.)-1;}
	    else{return(0);
	    }
    }
  if(position.y < (ratio_y*DCA+sp->DCU_length))
    {*y = (int)(position.y/sp->ypixelwidth +1.)-1;}
     else{if(position.y > (ratio_y*DCA+sp->DCU_length+sp->DCU_gap)
			   && position.y < (ratio_y*DCA+DCA_length))
	 {*y = (int)(position.y/sp->ypixelwidth +1.)-1;}
	    else{return(0);
	    }
    }
  }else{//detector without gaps
    *x = (int)(position.x/sp->xpixelwidth +1.)-1;
    *y = (int)(position.y/sp->ypixelwidth +1.)-1;
  }

  if ((*x>=0) && (*x<sp->xwidth) && (*y>=0) && (*y<sp->ywidth)) {
    return (1); // Valid pixel.
  } else {
    return (0); // Invalid pixel.
  }
}


int getSquarePixel_protoMirax(SquarePixels* sp, struct Point2d position, int* x, int* y)
{
  int ratio_x=(int)((position.x)/sp->DCU_length);
  int ratio_y=(int)((position.y)/sp->DCU_length);

  if(position.x > (ratio_x*sp->DCU_length+sp->DCU_gap) && position.x < (ratio_x*sp->DCU_length+sp->DCU_gap+2*sp->xpixelwidth))
    {*x = (int)(position.x/sp->xpixelwidth +1.)-1;}
  else{
    return(0);
  }
  if(position.y > (ratio_y*sp->DCU_length+sp->DCU_gap) && position.y < (ratio_y*sp->DCU_length+sp->DCU_gap+2*sp->ypixelwidth))
    {*y = (int)(position.y/sp->ypixelwidth +1.)-1;}
  else{
    return(0);
  }

  if ((*x>=0) && (*x<sp->xwidth) && (*y>=0) && (*y<sp->ywidth)) {
    return (1); // Valid pixel.
  } else {
    return (0); // Invalid pixel.
  }
}





void SPupdateValidFlag(SquarePixels* sp, int* x, int* y, int nsplits)
{
  int valid_flag = 1;

  // First check the affected pixels and the surrounding pixels
  // whether there are already any events.
  int split;
  int xcount, ycount;
  for (split=0; split<nsplits; split++) {
    if (x[split] != INVALID_PIXEL) {
      for (xcount=MAX(x[split]-1,0); xcount<MIN(x[split]+2,sp->xwidth); xcount++) {
	for (ycount=MAX(y[split]-1,0); ycount<MIN(y[split]+2,sp->ywidth); ycount++) {
	  if (0 != sp->array[xcount][ycount].valid_flag) {
	    valid_flag = -1;
	    if (sp->array[xcount][ycount].valid_flag > 0) {
	      SPsetInvalidFlag(sp, xcount, ycount);
	    }
	  }
	}
      }
    }
  }

  // Update the affected pixels with the appropriate valid flag.
  for (split=0; split<nsplits; split++) {
    if (x[split] != INVALID_PIXEL) {
      sp->array[x[split]][y[split]].valid_flag = valid_flag;
    }
  }
}



void SPsetInvalidFlag(SquarePixels* sp, int x, int y)
{
  sp->array[x][y].valid_flag = -1;

  int xcount, ycount;
  for (xcount=MAX(x-1,0); xcount<MIN(x+2,sp->xwidth); xcount++) {
    for (ycount=MAX(y-1,0); ycount<MIN(y+2,sp->ywidth); ycount++) {
      if (sp->array[xcount][ycount].valid_flag > 0) {
	SPsetInvalidFlag(sp, xcount, ycount);
      }
    }
  }
}



void SPaddCharge(SquarePixels* sp, int x, int y, float charge) {
  if ((x >= 0) && (x < sp->xwidth) && (y >= 0) && (y < sp->ywidth)) {
    // Add the charge [keV] to the pixel.
    sp->array[x][y].charge += charge;

    // Set the read-out flag for this line. That means that the line has
    // been affected by an event since the last read-out cycle. Some kinds
    // of detector models require this information to determine, whether
    // it is necessary to read out the line in the next cycle.
    sp->line2readout[y] = 1;
  }
}
