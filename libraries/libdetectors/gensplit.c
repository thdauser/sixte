#include "gensplit.h"


////////////////////////////////////////////////////////////////////
// Static function declarations
////////////////////////////////////////////////////////////////////


/** Determines the minimum distance value out of an array with 4
    entries and returns the corresponding index. */
static inline int getMinimumDistance(double array[]);


////////////////////////////////////////////////////////////////////
// Program Code
////////////////////////////////////////////////////////////////////


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
  // TODO Introduce pile-up flag.

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
      
    // Add the charge (photon energy) to the affected pixel.
    addGenDetCharge2Pixel(detline[line], column, charge);

    
  } else if (GS_GAUSS==split->type) {  
    // Gaussian split model.

    // Number of affected pixels.
    int npixels = 0; 
    // x- and y-indices of affected pixels.
    int x[4], y[4];
    // Charge fractions in the individual pixels.
    float fraction[4];
    // Charge cloud size (3 sigma).
    const float ccsize = split->par1*3.;

    // The following array entries are used to transform between 
    // different array indices.
    int xe[4] = {1, 0,-1, 0};   
    int ye[4] = {0, 1, 0,-1};

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

    // Add charge to all valid pixels of the split event.
    int count;
    for(count=0; count<npixels; count++) {
      if ((x[count]>=0) || (y[count]>=0)) {
	addGenDetCharge2Pixel(detline[y[count]], x[count], charge*fraction[count]);
      }
    }

    // END of Gaussian split model.

  } else {
    printf("Error: split model not supported!");
    exit(0);
  }


  // TODO Call the event trigger routine.

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


