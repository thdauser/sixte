#include "phproj.h"


void phproj(GenDet* const det,
	    AttitudeCatalog* const ac,
	    PatternFile* const plf,
	    const double t0,
	    const double exposure,
	    int* const status)
{
  // LOOP over all patterns in the FITS table.
  long row;
  for (row=0; row<plf->nrows; row++) {
    
    // Read the next pattern from the file.
    Pattern pattern;
    getPatternFromFile(plf, row+1, &pattern, status);
    CHECK_STATUS_BREAK(*status);

    // Check whether we are still within the requested time interval.
    if (pattern.time < t0) continue;
    if (pattern.time > t0+exposure) break;

    // Determine the Position of the source on the sky:
    // First determine telescope pointing direction at the current time.
    Vector nx, ny, nz;
    getTelescopeAxes(ac, &nx, &ny, &nz, pattern.time, status);
    CHECK_STATUS_BREAK(*status);

    // Determine RA and DEC of the photon origin.
    // Exact position on the detector:
    struct Point2d detpos;
    detpos.x = // in [m]
      (pattern.rawx*1.-det->pixgrid->xrpix+0.5+sixt_get_random_number())*
      det->pixgrid->xdelt + 
      det->pixgrid->xrval;
    detpos.y = // in [m]
      (pattern.rawy*1.-det->pixgrid->yrpix+0.5+sixt_get_random_number())*
      det->pixgrid->ydelt + 
      det->pixgrid->yrval;
    
    // Determine the source position on the sky using the telescope 
    // axis pointing vector and a vector from the point of the intersection 
    // of the optical axis with the sky plane to the source position.
    Vector srcpos;
    srcpos.x = nz.x 
      +detpos.x/det->focal_length*nx.x
      +detpos.y/det->focal_length*ny.x;
    srcpos.y = nz.y 
      +detpos.x/det->focal_length*nx.y
      +detpos.y/det->focal_length*ny.y;
    srcpos.z = nz.z 
      +detpos.x/det->focal_length*nx.z
      +detpos.y/det->focal_length*ny.z;
    srcpos = normalize_vector(srcpos);

    // Determine the equatorial coordinates RA and DEC
    // (RA and DEC are in the range [-pi:pi] and [-pi/2:pi/2] respectively).
    calculate_ra_dec(srcpos, &pattern.ra, &pattern.dec);
    
    // Update the data in the pattern list FITS file.
    updatePatternInFile(plf, row+1, &pattern, status);
    CHECK_STATUS_BREAK(*status);
  } 
  CHECK_STATUS_VOID(*status);
  // END of LOOP over all patterns.
}
