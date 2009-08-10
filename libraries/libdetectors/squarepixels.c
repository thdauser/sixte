#include "squarepixels.h"


int initSquarePixels(SquarePixels* sp, struct SquarePixelsParameters* parameters) 
{
  int status = EXIT_SUCCESS;

  // Set the array dimensions:
  sp->xwidth = parameters->xwidth;
  sp->ywidth = parameters->ywidth;
  sp->xoffset = parameters->xwidth/2;
  sp->yoffset = parameters->ywidth/2;
  sp->xpixelwidth = parameters->xpixelwidth;
  sp->ypixelwidth = parameters->ypixelwidth;

  // Get the memory for the pixels:
  sp->array = (SquarePixel**)malloc(sp->xwidth*sizeof(SquarePixel*));
  if (NULL==sp->array) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: memory allocation for pixel array failed!\n", status);
    return(status);
  } else {
    int count;
    for(count=0; count<sp->xwidth; count++) {
      sp->array[count] = (SquarePixel*)malloc(sp->ywidth*sizeof(SquarePixel));
      if (NULL==sp->array[count]) {
	status = EXIT_FAILURE;
	HD_ERROR_THROW("Error: memory allocation for pixel array failed!\n", status);
	return(status);
      }
    }
  }

  // Clear the pixels.
  clearSquarePixels(sp);

  return(status);
}



////////////////////////////////////////////////////////////////
inline void clearLineSquarePixels(SquarePixels* sp, const int line) 
{
  int x;
  for (x=0; x<sp->xwidth; x++) {
    sp->array[x][line].charge = 0.;
  }
}




////////////////////////////////////////////////////////////////
inline void clearSquarePixels(SquarePixels* sp) 
{
  int line;
  for (line=0; line<sp->ywidth; line++) {
    clearLineSquarePixels(sp, line);
  }
}




////////////////////////////////////////////////////////////////
void cleanupSquarePixels(SquarePixels* sp) 
{
  if (NULL!=sp->array) {
    int count;
    for (count=0; count<sp->xwidth; count++) {
      if (NULL!=sp->array[count]) {
	free(sp->array[count]);
      }
    }
    free(sp->array);
    sp->array=NULL;
  }
}




////////////////////////////////////////////////////////////////
/** Determines the minimum distance value out of an array with 4
 * entries and returns the corresponding index. */
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




///////////////////////////////////////
int getSquarePixelsSplits(SquarePixels* sp, GenericDetector* gd, struct Point2d position, 
			  int* x, int* y, double* fraction)
{
  //#ifdef GAUSSIAN_CHARGE_CLOUDS
  int npixels = 0; // Number of affected pixels.

  // The following array entries are used to transform between 
  // different array indices.
  int xe[4] = {1,0,-1,0};   
  int ye[4] = {0,1,0,-1};

  // Calculate pixel indices (integer) of central affected pixel:
  x[0] = (int)(position.x/sp->xpixelwidth + (double)(sp->xwidth/2) +1.)-1;
  y[0] = (int)(position.y/sp->ypixelwidth + (double)(sp->ywidth/2) +1.)-1;
  
  // If charge cloud size is 0, i.e. no splits are created.
  if (gd->ccsize < 1.e-20) {
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
    if (distances[mindist] < gd->ccsize) {
      // Not a single event!
      x[1] = x[0] + xe[mindist];
      y[1] = y[0] + ye[mindist];

      double mindistgauss = gaussint(distances[mindist]/gd->ccsigma);

      double minimum = distances[mindist];
      // search for the next to minimum distance to an edge
      distances[mindist] = -1.;
      int secmindist = getMinimumDistance(distances);
      distances[mindist] = minimum;

      if (distances[secmindist] < gd->ccsize) {
	// Quadruple!
	npixels = 4;

	x[2] = x[0] + xe[secmindist];
	y[2] = y[0] + ye[secmindist];
	x[3] = x[1] + xe[secmindist];
	y[3] = y[1] + ye[secmindist];

	// Calculate the different charge fractions in the 4 affected pixels.
	double secmindistgauss = gaussint(distances[secmindist]/gd->ccsigma);
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

  // Check whether all pixels lie inside the detector:
  int count;
  for(count=0; count<npixels; count++) {
    if ((x[count]<0) || (x[count]>=sp->xwidth) ||
	(y[count]<0) || (y[count]>=sp->ywidth)) {
      x[count] = INVALID_PIXEL;
      y[count] = INVALID_PIXEL;
    }
  }

  // Return the number of affected pixels, or 0, if the position lies outside
  // the detector array.
  return(npixels);

  //#else 
  /*
  // None-Gaussian Charge clouds:
  // Calculate the carge splitting according to a concept proposed by Konrad Dennerl.

  int npixels = 0;  // number of affected pixels

  // The following array entries are used to transform between 
  // different array indices.
  int xe[4] = {1,0,-1,0};   
  int ye[4] = {0,1,0,-1};

  // -- Determine the affected pixels:

  // Calculate pixel indices (integer) of central affected pixel:
  x[0] = (int)(position.x/detector->pixelwidth + (double)(detector->width/2) +1.)-1;
  y[0] = (int)(position.y/detector->pixelwidth + (double)(detector->width/2) +1.)-1;
  
  // Calculate the distances from the event to the borders of the 
  // surrounding pixel (in [m]):
  double distances[4] = { 
    // distance to right pixel edge
    (x[0]-detector->offset+1)*detector->pixelwidth - position.x,
    // distance to upper edge
    (y[0]-detector->offset+1)*detector->pixelwidth - position.y,
    // distance to left pixel edge
    position.x - (x[0]-detector->offset)*detector->pixelwidth,
    // distance to lower edge
    position.y - (y[0]-detector->offset)*detector->pixelwidth
  };

  int mindist = min_dist(distances, 4);
  x[1] = x[0] + xe[mindist];
  y[1] = y[0] + ye[mindist];

  double minimum = distances[mindist];
  // search for the next to minimum distance to an edge
  distances[mindist] = -1.;
  int secmindist = min_dist(distances, 4);
  distances[mindist] = minimum;

  x[2] = x[0] + xe[secmindist];
  y[2] = y[0] + ye[secmindist];
  x[3] = x[1] + xe[secmindist];
  y[3] = y[1] + ye[secmindist];

  // --- Now we know the affected pixels and can determine the charge fractions.

  fraction[0] = exp(-pow((sqrt(pow(detector->pixelwidth/2-distances[mindist],2.)+
			       pow(detector->pixelwidth/2-distances[secmindist],2.)))
			 /(0.35*75.e-6),2.));
  fraction[1] = exp(-pow((sqrt(pow(detector->pixelwidth/2+distances[mindist],2.)+
			       pow(detector->pixelwidth/2-distances[secmindist],2.)))
			 /(0.35*75.e-6),2.));
  fraction[2] = exp(-pow((sqrt(pow(detector->pixelwidth/2-distances[mindist],2.)+
			       pow(detector->pixelwidth/2+distances[secmindist],2.)))
			 /(0.35*75.e-6),2.));
  fraction[3] = exp(-pow((sqrt(pow(detector->pixelwidth/2+distances[mindist],2.)+
			       pow(detector->pixelwidth/2+distances[secmindist],2.)))
			 /(0.35*75.e-6),2.));
  double sum = fraction[0]+fraction[1]+fraction[2]+fraction[3];
  fraction[0] /= sum;
  fraction[1] /= sum;
  fraction[2] /= sum;
  fraction[3] /= sum;


  npixels = 4;


  // --- Check whether all pixels lie inside the detector:
  int count;
  for(count=0; count<npixels; count++) {
    if ((x[count]<0) || (x[count]>=detector->width) ||
	(y[count]<0) || (y[count]>=detector->width)) {
      x[count] = INVALID_PIXEL;
      y[count] = INVALID_PIXEL;
    }
  }

  // Return the number of affected pixel, and 0, if the position lies outside
  // the detector array.
  return(npixels);
  */
  //#endif /* GAUSSIAN_CHARGE_CLOUDS */
}

