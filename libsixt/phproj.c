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


void phproj_advdet(GenInst* const inst,
		AdvDet* const adv_det,
	    Attitude* const ac,
	    TesEventFile* const event_file,
	    const double t0,
	    const double exposure,
	    char proj_center,
	    int* const status)
{
	//const double cosrota=cos(inst->det->pixgrid->rota);
	//const double sinrota=sin(inst->det->pixgrid->rota);


	// Load detector geometry from advanced XML file
	Obj2D_instance * general_geometry= getObj2DFromAdvdet(adv_det,status);
	CHECK_STATUS_VOID(*status);

	// LOOP over all events in the input file.
	long pixid;
	double time,ra,dec;
	int anynul=0;
	while (event_file->row <= event_file->nrows) {

		// Read the next pixID from the file.
		fits_read_col(event_file->fptr, TLONG, event_file->pixIDCol,
							  event_file->row,1,1,0,&pixid, &anynul,status);
		fits_read_col(event_file->fptr, TDOUBLE, event_file->timeCol,
							  event_file->row,1,1,0,&time, &anynul,status);
		CHECK_STATUS_BREAK(*status);

		// Check whether we are still within the requested time interval.
		if (time < t0) continue;
		if (time > t0+exposure) break;

		// Determine the Position of the source on the sky:
		// First determine telescope pointing direction at the current time.
		Vector nx, ny, nz;
		getTelescopeAxes(ac, &nx, &ny, &nz, time, status);
		CHECK_STATUS_BREAK(*status);

		// Determine RA and DEC of the photon origin.
		// Exact position on the pixel.
		double xpix=0.,ypix=0.;
		if (!proj_center){
			xpix=(sixt_get_random_number(status)-0.5)*general_geometry->subobj[pixid-1]->geometry->width;
			CHECK_STATUS_BREAK(*status);
			ypix=(sixt_get_random_number(status)-0.5)*general_geometry->subobj[pixid-1]->geometry->height;
			CHECK_STATUS_BREAK(*status);
		}
		// Exact position in the detector sytem
//		double xr=cosrota*xb -sinrota*yb;
//		double yr=sinrota*xb +cosrota*yb;
		double cosrota = cos(general_geometry->subobj[pixid-1]->geometry->rota);
		double sinrota = sin(general_geometry->subobj[pixid-1]->geometry->rota);
		struct Point2d detpos;
		detpos.x=general_geometry->subobj[pixid-1]->geometry->cx + cosrota*xpix - sinrota*ypix;
		detpos.y=general_geometry->subobj[pixid-1]->geometry->cy + sinrota*xpix + cosrota*ypix;

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
		calculate_ra_dec(srcpos, &ra, &dec);

		// Update the data in the TES event list file.
		updateRaDecDetXY(event_file,ra,dec,detpos.x,detpos.y,status);
		CHECK_STATUS_BREAK(*status);

		event_file->row++;
	}
	CHECK_STATUS_VOID(*status);
	freeObj2D_instance(general_geometry);
	free(general_geometry);
	general_geometry=NULL;
	// END of LOOP over all events.
}
