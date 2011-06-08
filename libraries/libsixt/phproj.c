#include "phproj.h"


void phproj(GenDet* const det,
	    AttitudeCatalog* const ac,
	    EventListFile* const elf,
	    const double t0,
	    const double exposure,
	    int* const status)
{
  // LOOP over all events in the FITS table.
  long row;
  for (row=0; row<elf->nrows; row++) {
    
    // Read the next event from the file.
    Event event;
    getEventFromFile(elf, row+1, &event, status);
    CHECK_STATUS_BREAK(*status);

    // Check whether we are still within the requested time interval.
    if (event.time < t0) continue;
    if (event.time > t0+exposure) break;

    // Determine the Position of the source on the sky:
    // First determine telescope pointing direction at the current time.
    Vector nx, ny, nz;
    getTelescopeAxes(ac, &nx, &ny, &nz, event.time, status);
    CHECK_STATUS_BREAK(*status);

    // Determine RA and DEC of the photon origin.
    // Exact position on the detector:
    struct Point2d detpos;
    detpos.x = // in [m]
      (event.rawx*1.-det->pixgrid->xrpix+0.5+sixt_get_random_number())*
      det->pixgrid->xdelt + 
      det->pixgrid->xrval;
    detpos.y = // in [m]
      (event.rawy*1.-det->pixgrid->yrpix+0.5+sixt_get_random_number())*
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
    calculate_ra_dec(srcpos, &event.ra, &event.dec);
    
    // Update the data in the Event List FITS file.
    updateEventInFile(elf, row+1, &event, status);
    CHECK_STATUS_BREAK(*status);
  } 
  CHECK_STATUS_VOID(*status);
  // END of LOOP over all events.
}
