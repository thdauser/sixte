#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <limits.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#include "sixt.h"
#include "vector.h"
#include "check_fov.h"
#include "fits_ctlg.h"
#include "random.h"
#include "psf.h"
#include "photon.h"
#include "telescope.h"
#include "orbatt.h"


#define TOOLSUB event_projection_main
#include "headas_main.c"


#define REFXCUNI "deg"         // WCS physical unit of X axis 
#define REFXCRPX (12960000)    // WCS axis reference pixel
#define REFXCRVL (0.)          // [deg] WCS coord. at X axis ref. pixel
#define REFXCDLT (1.38889e-05) // [deg/pix] WCS X increment at ref. pixel (0.05"/pixel)

#define REFYCUNI "deg"         // WCS  physical unit of Y axis 
#define REFYCRPX (6480000)     // WCS axis reference pixel
#define REFYCRVL (0.)          // [deg] WCS coord. at Y axis ref. pixel
#define REFYCDLT (1.38889e-05) // [deg/pix] WCS Y increment at ref. pixel (0.05"/pixel)


/* Program parameters */
struct Parameters {
  char attitude_filename[FILENAME_LENGTH];  // filename of the attitude file
  char eventlist_filename[FILENAME_LENGTH]; // input: photon list

  double focal_length;
  double fov_diameter;
};


int event_projection_getpar(struct Parameters *parameters);



////////////////////////////////////
/** Main procedure. */
int event_projection_main() {
  struct Parameters parameters;   // Program parameters

  AttitudeCatalog* attitudecatalog=NULL;
  struct Eventlist_File* eventlistfile;

  struct Telescope telescope; // Telescope data (like FOV diameter or focal length)

  char msg[MAXMSG];           // error output buffer
  int status=EXIT_SUCCESS;    // error status


  // Register HEATOOL:
  set_toolname("event_projection");
  set_toolversion("0.01");


  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // Read the program parameters using PIL library.
    if ((status=event_projection_getpar(&parameters))) break;

    // Based on the parameters set up the program configuration.
    telescope.focal_length = parameters.focal_length;
    telescope.fov_diameter = parameters.fov_diameter;

    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);

    // Open the FITS file with the input event list:
    eventlistfile=open_EventlistFile(parameters.eventlist_filename, READWRITE, &status);
    if ((EXIT_SUCCESS!=status)||(NULL==eventlistfile)) break;
    
    // Write header keywords.
    fits_write_key(eventlistfile->fptr, TSTRING, "REFXCUNI", REFXCUNI, 
		   "WCS physical unit of X axis", &status);
    long lbuffer = REFXCRPX;
    fits_write_key(eventlistfile->fptr, TLONG, "REFXCRPX", &lbuffer, 
		   "WCS axis reference pixel", &status);
    double dbuffer = REFXCRVL;
    fits_write_key(eventlistfile->fptr, TDOUBLE, "REFXCRVL", &dbuffer,
		   "[deg] WCS coord. at X axis ref. pixel", &status);
    dbuffer = REFXCDLT;
    fits_write_key(eventlistfile->fptr, TDOUBLE, "REFXCDLT", &dbuffer,
		   "[deg/pix] WCS X increment at ref. pixel", &status);

    fits_write_key(eventlistfile->fptr, TSTRING, "REFYCUNI", REFYCUNI, 
		   "WCS physical unit of Y axis", &status);
    lbuffer = REFYCRPX;
    fits_write_key(eventlistfile->fptr, TLONG, "REFYCRPX", &lbuffer, 
		   "WCS axis reference pixel", &status);
    dbuffer = REFYCRVL;
    fits_write_key(eventlistfile->fptr, TDOUBLE, "REFYCRVL", &dbuffer,
		   "[deg] WCS coord. at Y axis ref. pixel", &status);
    dbuffer = REFYCDLT;
    fits_write_key(eventlistfile->fptr, TDOUBLE, "REFYCDLT", &dbuffer,
		   "[deg/pix] WCS Y increment at ref. pixel", &status);


    // Determine the time of the first and of the last event in the list. 
    // (This data is needed to read the adequate orbit/attitude information.)
    // Read the first event from the FITS file.
    int anynul = 0.;
    double t0=0., timespan=0.;
    fits_read_col(eventlistfile->fptr, TDOUBLE, 1, 1, 1, 1, &t0, &t0, &anynul, &status);
    fits_read_col(eventlistfile->fptr, TDOUBLE, 1, eventlistfile->nrows, 1, 1, 
		  &timespan, &timespan, &anynul, &status);
    

    // Get the satellite catalog with the orbit and (telescope) attitude data:
    /* OBSOLETE
    if ((status=get_satellite_catalog(&sat_catalog, &sat_nentries, 
				      0., timespan+100., 
				      parameters.orbit_filename, 
				      parameters.attitude_filename))
	!=EXIT_SUCCESS) break;
    */
    if (NULL==(attitudecatalog=get_AttitudeCatalog(parameters.attitude_filename,
						   t0, timespan, &status))) break;
						   

    
    // Create an array that contains the off-axis angle corresponding to 
    // some particular positions on the detector.
    struct angle_position_relation { 
      double distance; // [m]
      double angle;    // offaxis-angle [rad]
    };
    struct angle_position_relation apr[7] = {
      { 0.      ,  0.                  },
      { 0.002338,  5. /60.*M_PI/180. },
      { 0.004662, 10. /60.*M_PI/180. },
      { 0.006988, 15. /60.*M_PI/180. },
      { 0.009288, 20. /60.*M_PI/180. },
      { 0.011613, 25. /60.*M_PI/180. },
      { 0.013938, 30. /60.*M_PI/180. }
    };
      
    
    // --- END of Initialization ---


    // --- Beginning of the Imaging Process ---

    // Beginning of actual simulation (after loading required data):
    headas_chat(5, "start imaging process ...\n");

    // LOOP over all event in the FITS table
    long attitude_counter=0;   // counter for entries in the AttitudeCatalog

    // SCAN PHOTON LIST    
    for(eventlistfile->row=0; 
	(eventlistfile->row<eventlistfile->nrows)&&(status==EXIT_SUCCESS); 
	eventlistfile->row++) {
      
      // Read the event from the FITS file.
      struct Event event;
      if (get_eventlist_row(*eventlistfile, &event, &status)) break;

      // Get the last orbit entry before 'event.time'
      // (in order to interpolate the attitude at this time between 
      // the neighboring calculated values):
      for( ; attitude_counter<attitudecatalog->nentries-1; attitude_counter++) {
	if(attitudecatalog->entry[attitude_counter+1].time>event.time) {
	  break;
	}
      }
      if(fabs(attitudecatalog->entry[attitude_counter].time-event.time)>60.) { 
	// no entry within 10 minutes !!
	status = EXIT_FAILURE;
	sprintf(msg, "Error: no adequate orbit entry for time %lf!\n", event.time);
	HD_ERROR_THROW(msg,status);
	break;
      }

      // Determine the Position of the source on the sky:

      // First determine telescope pointing direction at the current time.
      telescope.nz = 
	normalize_vector(interpolate_vec(attitudecatalog->entry[attitude_counter].nz, 
					 attitudecatalog->entry[attitude_counter].time, 
					 attitudecatalog->entry[attitude_counter+1].nz, 
					 attitudecatalog->entry[attitude_counter+1].time, 
					 event.time));
      // Determine the current nx: perpendicular to telescope axis nz
      // and in the direction of the satellite motion.
      telescope.nx = 
	normalize_vector(interpolate_vec(attitudecatalog->entry[attitude_counter].nx, 
					 attitudecatalog->entry[attitude_counter].time, 
					 attitudecatalog->entry[attitude_counter+1].nx, 
					 attitudecatalog->entry[attitude_counter+1].time, 
					 event.time));
      // Remove the component along the vertical direction nz 
      // (nx must be perpendicular to nz!):
      double scp = scalar_product(&telescope.nz, &telescope.nx);
      telescope.nx.x -= scp*telescope.nz.x;
      telescope.nx.y -= scp*telescope.nz.y;
      telescope.nx.z -= scp*telescope.nz.z;
      telescope.nx = normalize_vector(telescope.nx);
      // The third axis of the coordinate system ny is perpendicular 
      // to telescope axis nz and nx:
      telescope.ny=normalize_vector(vector_product(telescope.nz, telescope.nx));



      // Determine RA, DEC and the sky coordinates (in pixel) of the source.
      struct Point2d detector_position;
      detector_position.x = ((double)(event.xi-384/2)+0.5)*75.e-6; // in [m]
      detector_position.y = ((double)(event.yi-384/2)+0.5)*75.e-6; // in [m]
      double d = sqrt(pow(detector_position.x,2.)+pow(detector_position.y,2.));

      // Determine the offaxis_angle corresponding to the detector position.
      double offaxis_angle; // = atan(d/telescope.focal_length);
      // Interpolation:
      int count;
      for(count=1; count<7; count++) {
	if (apr[count].distance>d) break;
      }
      if (count<6) {
	offaxis_angle = apr[count-1].angle + 
	  (d-apr[count-1].distance)/(apr[count].distance-apr[count-1].distance) * 
	  (apr[count].angle-apr[count-1].angle);
      } else {
	offaxis_angle = apr[6].angle + 
	  (d-apr[6].distance)/(apr[6].distance-apr[5].distance) *
	  (apr[6].angle-apr[5].angle);
      }

      double r = tan(offaxis_angle);

      struct vector source_position;
      source_position.x = telescope.nz.x 
	-r*(detector_position.x/d*telescope.nx.x+detector_position.y/d*telescope.ny.x);
      source_position.y = telescope.nz.y 
	-r*(detector_position.x/d*telescope.nx.y+detector_position.y/d*telescope.ny.y);
      source_position.z = telescope.nz.z 
	-r*(detector_position.x/d*telescope.nx.z+detector_position.y/d*telescope.ny.z);
      source_position = normalize_vector(source_position);

      event.ra = atan2(source_position.y, source_position.x) * 180./M_PI;
      event.dec = asin(source_position.z) * 180./M_PI;

      // Put some randomization on the RA and DEC coordinate (within the sky pixel)
      // to receive a continuous image.
      // (spread by the width of one detector pixel on the sky)
      event.ra  += (get_random_number()-0.5)*0.00265625; // (= (61.2/60.)°/384)
      event.dec += (get_random_number()-0.5)*0.00265625;

      // Determine the pixel coordinates in the sky image:
      event.sky_xi = (int)((event.ra -REFXCRVL)/REFXCDLT+REFXCRPX);
      event.sky_yi = (int)((event.dec-REFYCRVL)/REFYCDLT+REFYCRPX);

      // Store the data in the Event List FITS file.
      fits_write_col(eventlistfile->fptr, TDOUBLE, 10, eventlistfile->row+1, 
		     1, 1, &event.ra, &status);
      fits_write_col(eventlistfile->fptr, TDOUBLE, 11, eventlistfile->row+1, 
		     1, 1, &event.dec, &status);
      fits_write_col(eventlistfile->fptr, TLONG, 12, eventlistfile->row+1, 
		     1, 1, &event.sky_xi, &status);
      fits_write_col(eventlistfile->fptr, TLONG, 13, eventlistfile->row+1, 
		     1, 1, &event.sky_yi, &status);


    } // END of scanning-LOOP over the event list.
  } while(0);  // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // release HEADAS random number generator
  HDmtFree();

  // Close the FITS files.
  if (eventlistfile->fptr) fits_close_file(eventlistfile->fptr, &status);

  // Release memory of orbit/attitude catalog
  //  if (sat_catalog) free(sat_catalog);

  free_AttitudeCatalog(attitudecatalog);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully!\n\n");

  return(status);
}




////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int event_projection_getpar(struct Parameters *parameters)
{
  char msg[MAXMSG];             // error output buffer
  int status=EXIT_SUCCESS;      // error status

  // Get the filename of the attitude file (FITS file)
  if ((status = PILGetFname("attitude_filename", parameters->attitude_filename))) {
    sprintf(msg, "Error reading the filename of the attitude file!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the filename of the input photon list (FITS file)
  else if ((status = PILGetFname("eventlist_filename", parameters->eventlist_filename))) {
    sprintf(msg, "Error reading the filename of the event list!\n");
    HD_ERROR_THROW(msg,status);
  }
  
  // Read the diameter of the FOV (in arcmin)
  else if ((status = PILGetReal("fov_diameter", &parameters->fov_diameter))) {
    sprintf(msg, "Error reading the diameter of the FOV!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Read the focal length [m]
  else if ((status = PILGetReal("focal_length", &parameters->focal_length))) {
    sprintf(msg, "Error reading the focal length!\n");
    HD_ERROR_THROW(msg,status);
  }

  // convert angles from [arc min] to [rad]
  parameters->fov_diameter = parameters->fov_diameter*M_PI/(60.*180.); 

  return(status);
}



