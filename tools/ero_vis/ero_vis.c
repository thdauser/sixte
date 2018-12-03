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

#include "sixt.h"
#include "vector.h"
#include "telescope.h"
#include "attitude.h"
#include "check_fov.h"
#include "gti.h"
#include "simput.h"

#define TOOLSUB ero_vis_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  /** Attitude file. */
  char* Attitude;

  /** Source catalog. */
  char* Simput;
  double SrcRA, SrcDec;

  /** GTI file. */
  char* GTIfile;

  double TSTART;
  double Exposure;
  double dt;

  /** [rad]. */
  double visibility_range;

  int clobber;
};


int ero_vis_getpar(struct Parameters *parameters);


int ero_vis_main()
{
  // Program parameters.
  struct Parameters par;

  Attitude* ac=NULL;
  SimputCtlg* cat=NULL;
  GTI* gti=NULL;
  fitsfile* fptr=NULL;

  char* datestr=NULL;
  char* timestr=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL:
  set_toolname("ero_vis");
  set_toolversion("0.08");


  do { // Beginning of the ERROR handling loop.

    // --- Initialization ---

    // Read the program parameters using PIL library.
    if ((status=ero_vis_getpar(&par))) break;

    // Get the telescope attitude data.
    ac=loadAttitude(par.Attitude, &status);
    CHECK_STATUS_BREAK(status);

    // Load the SIMPUT source catalog, if available.
    if ( par.Simput!=NULL ) {
    	cat=openSimputCtlg(par.Simput, READONLY, 0, 0, 0, 0, &status);
    	CHECK_STATUS_BREAK(status);

    	// Check if the catalog contains any sources.
    	if (0==cat->nentries) {
    		SIXT_ERROR("SIMPUT catalog is empty");
    		status=EXIT_FAILURE;
    		break;
    	}
    }
    // Otherwise use the specified RA and Dec source position.

    // Set up a new GTI collection.
    gti=newGTI(&status);
    CHECK_STATUS_BREAK(status);
    gti->mjdref=ac->mjdref;

    // --- END of Initialization ---


    // --- Beginning of source localization calculation ---

    headas_chat(3, "determine the search cone ...\n");

    // Determine a reference vector and a cone opening angle
    // including all sources in the catalog.
    Vector refpos;
    double cone_radius=0.; // [rad]

    // If a SIMPUT catalog has been specified.
    if (NULL!=cat) {
    	// Get the first source in the catalog.
    	SimputSrc* src=getSimputSrc(cat, 1, &status);
    	CHECK_STATUS_BREAK(status);

    	// Determine its position and angular extension.
    	refpos=unit_vector(src->ra, src->dec);
    	cone_radius=getSimputSrcExt(cat, src, 0., 0., &status);
    	CHECK_STATUS_BREAK(status);

    	// Loop over all sources in the catalog.
    	long ii;
    	for (ii=1; ii<cat->nentries; ii++) {
    		// Get the next source in the catalog.
    		SimputSrc* src=getSimputSrc(cat, ii+1, &status);
    		CHECK_STATUS_BREAK(status);

    		// Determine its position and angular extension.
    		Vector srcpos=unit_vector(src->ra, src->dec);
    		float extension=getSimputSrcExt(cat, src, 0., 0., &status);
    		CHECK_STATUS_BREAK(status);

    		// Determine the angle between the reference direction and
    		// the source position.
    		double angular_distance=acos(scalar_product(&refpos, &srcpos));

    		// Check if the search radius has to be enlarged.
    		if (angular_distance+extension>cone_radius) {
    			double delta1=extension+angular_distance-cone_radius;

    			// If the difference is small, simply enlarge the opening
    			// angle of the search cone.
    			if (delta1 < 1./3600.*M_PI/180.) {
    				cone_radius+=delta1;
    			} else {
    				double delta2=extension-angular_distance-cone_radius;
    				if (delta2<0.) {
    					delta2=0.;
    				}
    				refpos=interpolateCircleVector(refpos, srcpos,
    						(delta1-delta2)*0.5/angular_distance);
    				cone_radius+=0.5*(delta1+delta2);
    			}
    		}
    	}
    } else {
    	// Determine the reference position from the explicitly
    	// specified RA and Dec.
    	refpos=unit_vector(par.SrcRA *M_PI/180., par.SrcDec*M_PI/180.);
    }


    // Print some informational data.
    double ra, dec;
    calculate_ra_dec(refpos, &ra, &dec);
    headas_chat(5, "ra=%lf deg, dec=%lf deg, cone radius=%lf deg\n",
    		ra*180./M_PI, dec*180./M_PI, cone_radius*180./M_PI);

    // Determine the diameter of the search radius (minimum cos-value).
    // (angle(telescope,source) <= 1/2 * diameter)
    double search_angle=0.5*par.visibility_range+cone_radius;
    double min_align;
    if (search_angle<=M_PI) {
    	min_align=cos(search_angle);
    } else {
    	min_align=-1.;
    }

    // --- Beginning of GTI calculation ---

    headas_chat(3, "calculate the visibility GTIs ...\n");

    // LOOP over the given time interval in steps of dt.
    double time;
    double start=0;
    double ininterval=0;
    for (time=par.TSTART; time<par.TSTART+par.Exposure; time+=par.dt) {

    	// Print the current time (program status information for the user).
    	headas_chat(5, "\rtime: %.1lf s ", time);
    	fflush(NULL);

    	// Determine the telescope pointing direction at the current time.
    	Vector telescope_nz=getTelescopeNz(ac, time, &status);
    	CHECK_STATUS_BREAK(status);

    	// Check if the FOV touches the search cone.
    	if (check_fov(&refpos, &telescope_nz, min_align)==0) {
    		// Source lies inside the FOV.
    		if (0==ininterval) {
    			ininterval=1;
    			start=time;
    		}
    	} else {
    		// Source lies outside the FOV.
    		if (1==ininterval) {
    			ininterval=0;
    			appendGTI(gti, start, time, &status);
    			CHECK_STATUS_BREAK(status);
    		}
    	}
    }
    CHECK_STATUS_BREAK(status);
    // END of LOOP over the specified time interval.


    // Store the GTI in the output file.
    saveGTI(gti, par.GTIfile, par.clobber, &status);
    CHECK_STATUS_BREAK(status);


    // Open the GTI file and append the columns 'DATE-START', 'TIME-START',
    // 'DATE-STOP', and 'TIME-STOP', which contain the same information
    // as the already present columns 'START' and 'STOP', but are better
    // readable for human being than large numbers of seconds.

    // Open the file.
    fits_open_file(&fptr, par.GTIfile, READWRITE, &status);
    CHECK_STATUS_BREAK(status);
    int hdutype;
    fits_movabs_hdu(fptr, 2, &hdutype, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the number of rows.
    long nrows;
    fits_get_num_rows(fptr, &nrows, &status);
    CHECK_STATUS_BREAK(status);

    // Insert the new columns.
    fits_insert_col(fptr, 3, "DATE-START", "10A", &status);
    fits_insert_col(fptr, 4, "TIME-START", "8A", &status);
    fits_insert_col(fptr, 5, "DATE-STOP", "10A", &status);
    fits_insert_col(fptr, 6, "TIME-STOP", "8A", &status);
    CHECK_STATUS_BREAK(status);

    // Loop over all entries.
    datestr=(char*)malloc(20*sizeof(char));
    CHECK_NULL_BREAK(datestr, status, "memory allocation for string buffer failed");
    timestr=(char*)malloc(20*sizeof(char));
    CHECK_NULL_BREAK(datestr, status, "memory allocation for string buffer failed");
    long jj;
    for (jj=0; jj<nrows; jj++) {
      // Determine the start date and time.
      sixt_get_date_time(gti->mjdref, gti->start[jj], datestr, timestr, &status);
      CHECK_STATUS_BREAK(status);
      fits_write_col(fptr, TSTRING, 3, jj+1, 1, 1, &datestr, &status);
      fits_write_col(fptr, TSTRING, 4, jj+1, 1, 1, &timestr, &status);
      CHECK_STATUS_BREAK(status);

      // Determine the stop date and time.
      sixt_get_date_time(gti->mjdref, gti->stop[jj], datestr, timestr, &status);
      CHECK_STATUS_BREAK(status);
      fits_write_col(fptr, TSTRING, 5, jj+1, 1, 1, &datestr, &status);
      fits_write_col(fptr, TSTRING, 6, jj+1, 1, 1, &timestr, &status);
      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);

    // Close the file.
    fits_close_file(fptr, &status);
    CHECK_STATUS_BREAK(status);
    fptr=NULL;

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Close the GTI file.
  if (NULL!=fptr) {
    fits_close_file(fptr, &status);
  }

  // Release memory.
  if (NULL!=datestr) free(datestr);
  if (NULL!=timestr) free(timestr);
  freeAttitude(&ac);
  freeGTI(&gti);

  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully!\n\n");
  return(status);
}


int ero_vis_getpar(struct Parameters *par)
{
  int status=EXIT_SUCCESS; // Error status

  query_simput_parameter_file_name("Attitude", &par->Attitude, &status);
  query_simput_parameter_file_name("Simput", &par->Simput, &status);

  // only load srcRA and srcDec if Simput is not given
  if ( par->Simput==NULL ) {
	  query_simput_parameter_double("SrcRA",&par->SrcRA,&status);
	  query_simput_parameter_double("SrcDec",&par->SrcDec,&status);
	  headas_chat(3, "using SrcRA=%.3f, SrcDec=%.3f as no Simput file is given\n",par->SrcRA,par->SrcDec);
  } else {
	  // set to default values
	  par->SrcRA=0.0;
	  par->SrcDec=0.0;
	  headas_chat(3, "using Simput File: %s \n",par->Simput);
  }

  query_simput_parameter_file_name("GTIfile", &par->GTIfile, &status);
  query_simput_parameter_double("visibility_range", &par->visibility_range, &status);
  // Convert FOV diameter from [deg] to [rad].
  par->visibility_range*=M_PI/180.;
  query_simput_parameter_double("TSTART", &par->TSTART, &status);
  query_simput_parameter_double("Exposure", &par->Exposure, &status);
  query_simput_parameter_double("dt", &par->dt, &status);
  query_simput_parameter_bool("clobber", &par->clobber, &status);

  return(status);
}
