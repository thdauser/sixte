#include "phproj.h"


void phproj(GenInst* const inst,
	    Attitude* const ac,
	    EventFile* const elf,
	    const double t0,
	    const double exposure,
	    int* const status)
{
  const double cosrota=cos(inst->det->pixgrid->rota);
  const double sinrota=sin(inst->det->pixgrid->rota);

  // LOOP over all events in the input file.
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
    // Exact position on the detector.
    double xb=
      (event.rawx*1.-inst->det->pixgrid->xrpix
       +0.5+sixt_get_random_number(status))*inst->det->pixgrid->xdelt;
    CHECK_STATUS_BREAK(*status);
    double yb=
      (event.rawy*1.-inst->det->pixgrid->yrpix
       +0.5+sixt_get_random_number(status))*inst->det->pixgrid->ydelt;
    CHECK_STATUS_BREAK(*status);
    
    double xr=cosrota*xb -sinrota*yb;
    double yr=sinrota*xb +cosrota*yb;
    
    struct Point2d detpos;
    detpos.x=xr + inst->det->pixgrid->xrval; // in [m]      
    detpos.y=yr + inst->det->pixgrid->yrval; // in [m]
    
    // Determine the source position on the sky using the telescope 
    // axis pointing vector and a vector from the point of the intersection 
    // of the optical axis with the sky plane to the source position.
    Vector srcpos;
    srcpos.x = nz.x 
      +detpos.x/inst->tel->focal_length*nx.x
      +detpos.y/inst->tel->focal_length*ny.x;
    srcpos.y = nz.y 
      +detpos.x/inst->tel->focal_length*nx.y
      +detpos.y/inst->tel->focal_length*ny.y;
    srcpos.z = nz.z 
      +detpos.x/inst->tel->focal_length*nx.z
      +detpos.y/inst->tel->focal_length*ny.z;
    srcpos = normalize_vector(srcpos);

    // Determine the equatorial coordinates RA and DEC
    // (RA and DEC are in the range [-pi:pi] and [-pi/2:pi/2] respectively).
    calculate_ra_dec(srcpos, &event.ra, &event.dec);
    
    // Update the data in the event list FITS file.
    updateEventInFile(elf, row+1, &event, status);
    CHECK_STATUS_BREAK(*status);
  } 
  CHECK_STATUS_VOID(*status);
  // END of LOOP over all events.
}
