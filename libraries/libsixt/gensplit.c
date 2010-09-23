#include "gensplit.h"


////////////////////////////////////////////////////////////////////
// Static function declarations
////////////////////////////////////////////////////////////////////


/** Determines the minimum distance value out of an array with 4
    entries and returns the corresponding index. */
static inline int getMinimumDistance(const double array[]);

/** Set the pile-up flag of the specified pixel and all surrounding
    pixels containing a positive charge by recursive function
    calls. */
static void setGenPileupFlag(GenDetLine** const line,	
			     const GenPixGrid* const grid,
			     const int x, const int y);


////////////////////////////////////////////////////////////////////
// Program Code
////////////////////////////////////////////////////////////////////


GenSplit* newGenSplit(int* const status) 
{
  // Allocate memory.
  GenSplit* split=(GenSplit*)malloc(sizeof(GenSplit));
  if (NULL==split) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenSplit failed!\n", 
		   *status);
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
  // Number of affected pixels.
  int npixels=0;
  // x- and y-indices of affected pixels.
  int x[4], y[4];
  // Charge fractions in the individual pixels.
  float fraction[4];

  // The following array entries are used to transform between 
  // different array indices for accessing neighboring pixels.
  const int xe[4] = {1, 0,-1, 0};
  const int ye[4] = {0, 1, 0,-1};

  // Which kind of split model has been selected?
  if (GS_NONE==split->type) {
    // No split events => all events are singles.
    npixels=1;

    // Determine the affected detector line and column.
    x[0] = getGenDetAffectedColumn(grid, position->x);
    y[0] = getGenDetAffectedLine  (grid, position->y);

    // Check if the returned values are valid line and column indices.
    if ((x[0]<0) || (y[0]<0)) {
      return;
    }
    
  } else if (GS_GAUSS==split->type) {  
    // Gaussian split model.

    // Charge cloud size (3 sigma).
    const float ccsize = split->par1*3.;

    // Calculate pixel indices (integer) of the central affected pixel:
    x[0] = getGenDetAffectedColumn(grid, position->x);
    y[0] = getGenDetAffectedLine  (grid, position->y);
  
    // Check if the impact position lies inside the detector pixel array.
    if ((0>x[0]) || (0>y[0])) {
      return;
    }

    // Calculate the distances from the impact center position to the 
    // borders of the surrounding pixel (in [m]):
    double distances[4] = {
      // Distance to right pixel edge.
      (x[0]-grid->xrpix+1.5)*grid->xdelt + grid->xrval - position->x,
      // Distance to upper edge.
      (y[0]-grid->yrpix+1.5)*grid->ydelt + grid->yrval - position->y,
      // Distance to left pixel edge.
      position->x - ((x[0]-grid->xrpix+0.5)*grid->xdelt + grid->xrval),
      // distance to lower edge
      position->y - ((y[0]-grid->yrpix+0.5)*grid->ydelt + grid->yrval)
    };

    int mindist = getMinimumDistance(distances);
    if (distances[mindist] < ccsize) {
      // Not a single event!
      x[1] = x[0] + xe[mindist];
      y[1] = y[0] + ye[mindist];

      double mindistgauss = gaussint(distances[mindist]/split->par1);

      // Search for the next to minimum distance to an edge.
      double minimum = distances[mindist];
      distances[mindist] = -1.;
      int secmindist = getMinimumDistance(distances);
      distances[mindist] = minimum;

      if (distances[secmindist] < ccsize) {
	// Quadruple!
	npixels = 4;

	x[2] = x[0] + xe[secmindist];
	y[2] = y[0] + ye[secmindist];
	x[3] = x[1] + xe[secmindist];
	y[3] = y[1] + ye[secmindist];

	// Calculate the different charge fractions in the 4 affected pixels.
	double secmindistgauss = gaussint(distances[secmindist]/split->par1);
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
      
    } 
    // END of check for single event

    // END of Gaussian split model.

  } else if (GS_EXPONENTIAL==split->type) {  
    // Exponential split model.
    // None-Gaussian, exponential charge cloud model 
    // (concept proposed by Konrad Dennerl).
    npixels=4;

    // Calculate pixel indices (integer) of central affected pixel:
    x[0] = getGenDetAffectedColumn(grid, position->x);
    y[0] = getGenDetAffectedLine  (grid, position->y);
  
    // Check if the impact position lies inside the detector pixel array.
    if ((0>x[0]) || (0>y[0])) {
      return;
    }

    // Calculate the distances from the impact center position to the 
    // borders of the surrounding pixel (in units [fraction of a pixel edge]):
    double distances[4] = {
      // Distance to right pixel edge.
      x[0]-grid->xrpix+1.5 + (grid->xrval - position->x)/grid->xdelt,
      // Distance to upper edge.
      y[0]-grid->yrpix+1.5 + (grid->yrval - position->y)/grid->ydelt,
      // Distance to left pixel edge.
      (position->x - grid->xrval)/grid->xdelt - (x[0]-grid->xrpix+0.5),
      // distance to lower edge
      (position->y - grid->yrval)/grid->ydelt - (y[0]-grid->yrpix+0.5)
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
		      pow(split->par1,2.));
    fraction[1] = exp(-(pow(0.5+distances[mindist],2.)+pow(0.5-distances[secmindist],2.))/
		      pow(split->par1,2.));
    fraction[2] = exp(-(pow(0.5-distances[mindist],2.)+pow(0.5+distances[secmindist],2.))/
		      pow(split->par1,2.));
    fraction[3] = exp(-(pow(0.5+distances[mindist],2.)+pow(0.5+distances[secmindist],2.))/
		      pow(split->par1,2.));
    // Normalization to 1.
    double sum = fraction[0]+fraction[1]+fraction[2]+fraction[3];
    fraction[0] /= sum;
    fraction[1] /= sum;
    fraction[2] /= sum;
    fraction[3] /= sum;

    // END of exponential split model.

  } else {
    printf("Error: split model not supported!\n");
    exit(0);
  }


  // Manage the pile-up flags for the new charges.
  // Search for surrounding pixels already containing charges.
  int xmin = x[0];
  int ymin = y[0];
  int xmax = x[0];
  int ymax = y[0];
  int ii, jj;
  GenPileupFlag pileup = GP_NONE;
  for(ii=1; ii<npixels; ii++) {
    xmin = MIN(xmin, x[ii]);
    ymin = MIN(ymin, y[ii]);
    xmax = MAX(xmax, x[ii]);
    ymax = MAX(ymax, y[ii]);
  }
  xmin = MAX(0, xmin-1);
  ymin = MAX(0, ymin-1);
  xmax = MIN(grid->xwidth-1, xmax+1);
  ymax = MIN(grid->ywidth-1, ymax+1);
  for (ii=xmin; ii<=xmax; ii++) {
    for (jj=ymin; jj<=ymax; jj++) {
      if (detline[jj]->charge[ii]>0.) {
	pileup = GP_PILEUP;
      }
    }
  }
  // END of loop over all surrounding pixels.


  // Add charge to all valid pixels of the split event.
  for(ii=0; ii<npixels; ii++) {
    if ((x[ii]>=0) && (x[ii]<grid->xwidth) &&
	(y[ii]>=0) && (y[ii]<grid->ywidth)) {
      addGenDetCharge2Pixel(detline[y[ii]], x[ii], charge*fraction[ii]);
    }
  }


  // If necessary set the pile-up flag for the generated pattern.
  if (GP_PILEUP==pileup) {
    setGenPileupFlag(detline, grid, x[0], y[0]);
  }
  

  // TODO Call the event trigger routine.

}



static inline int getMinimumDistance(const double array[]) 
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



static void setGenPileupFlag(GenDetLine** const line,	
			     const GenPixGrid* const grid,
			     const int x, const int y)
{
  // Set the GenPileupFlag for the specified pixel.
  line[y]->pileup[x] = GP_PILEUP;

  // Check the surrounding pixels.
  int xmin = MAX(0, x-1);
  int ymin = MAX(0, y-1);
  int xmax = MIN(grid->xwidth-1, x+1);
  int ymax = MIN(grid->ywidth-1, y+1);
  int ii, jj;
  for (ii=xmin; ii<=xmax; ii++) {
    for (jj=ymin; jj<=ymax; jj++) {
      if (line[jj]->charge[ii] > 0.) {
	if (GP_PILEUP!=line[jj]->pileup[ii]) {
	  setGenPileupFlag(line, grid, ii, jj);
	}
      }
    }
  }
}

