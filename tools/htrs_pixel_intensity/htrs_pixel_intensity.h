#ifndef HTRS_PIXEL_INTENSITY_H
#define HTRS_PIXEL_INTENSITY_H 1

#include "sixt.h"
#include "hexagonalpixels.h"
#include "arcpixels.h"
#include "impact.h"
#include "impactlistfile.h"


#define TOOLSUB htrs_pixel_intensity_main
#include "headas_main.c"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////

struct Parameters {

  char impactlist_filename[FILENAME_LENGTH];

  /** Width of the spokes of the HTRS mask. */
  double mask_spoke_width;

};


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


// Reads the program parameters using PIL
static int getpar(struct Parameters* parameters);

#endif /* HTRS_PIXEL_INTENSITY_H */
