#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include "sixt.h"
#include "vector.h"
#include "point.h"
#include "sixt_random.h"
#include "telescope.h"
#include "attitudecatalog.h"
#include "erositaeventfile.h"

#define TOOLSUB event_projection_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char attitude_filename[FILENAME_LENGTH];  // filename of the attitude file
  char eventlist_filename[FILENAME_LENGTH]; // input: event list
  
  double t0;
  double timespan;

  double focal_length;
  double fov_diameter;
};


int event_projection_getpar(struct Parameters *parameters);


////////////////////////////////////
/** Main procedure. */
int event_projection_main() {
  struct Parameters parameters;   // Program parameters

  AttitudeCatalog* attitudecatalog=NULL;
  eROSITAEventFile eventlistfile;

  // WCS keywords:
  double tcrpxx, tcrvlx, tcdltx;
  double tcrpxy, tcrvly, tcdlty;

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
    telescope.fov_diameter = parameters.fov_diameter; // [rad]
    double radperpixel = telescope.fov_diameter/384;

    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);

    // Open the FITS file with the input event list:
    status=openeROSITAEventFile(&eventlistfile, parameters.eventlist_filename, READWRITE);
    if (EXIT_SUCCESS!=status) break;
    
    /*
    // Write header keywords.
    if (fits_write_key(eventlistfile.generic.fptr, TSTRING, "REFXCUNI", REFXCUNI, 
		       "WCS physical unit of X axis", &status)) break;
    long lbuffer = REFXCRPX;
    if (fits_write_key(eventlistfile.generic.fptr, TLONG, "REFXCRPX", &lbuffer, 
		       "WCS axis reference pixel", &status)) break;
    double dbuffer = REFXCRVL;
    if (fits_write_key(eventlistfile.generic.fptr, TDOUBLE, "REFXCRVL", &dbuffer,
		       "[deg] WCS coord. at X axis ref. pixel", &status)) break;
    dbuffer = REFXCDLT;
    if (fits_write_key(eventlistfile.generic.fptr, TDOUBLE, "REFXCDLT", &dbuffer,
		       "[deg/pix] WCS X increment at ref. pixel", &status)) break;

    if (fits_write_key(eventlistfile.generic.fptr, TSTRING, "REFYCUNI", REFYCUNI, 
		       "WCS physical unit of Y axis", &status)) break;
    lbuffer = REFYCRPX;
    if (fits_write_key(eventlistfile.generic.fptr, TLONG, "REFYCRPX", &lbuffer, 
		       "WCS axis reference pixel", &status)) break;
    dbuffer = REFYCRVL;
    if (fits_write_key(eventlistfile.generic.fptr, TDOUBLE, "REFYCRVL", &dbuffer,
		       "[deg] WCS coord. at Y axis ref. pixel", &status)) break;
    dbuffer = REFYCDLT;
    if (fits_write_key(eventlistfile.generic.fptr, TDOUBLE, "REFYCDLT", &dbuffer,
		       "[deg/pix] WCS Y increment at ref. pixel", &status)) break;
    */


    // Read HEADER keywords.
    char comment[MAXMSG]; // buffer
    // Attitude File:
    if (fits_read_key(eventlistfile.generic.fptr, TSTRING, "ATTITUDE", 
		      &parameters.attitude_filename, comment, &status)) break;
    if (0==strlen(parameters.attitude_filename)) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: no attitude file specified in FITS header of event list!\n",
		     status);
      break;
    }
    
    // WCS keywords:
    char columnnumber[3];
    char tcrpx[10];
    char tcrvl[10];
    char tcdlt[10];
    // For the sky x-coordinate:
    strcpy(tcrpx, "TCRPX");
    strcpy(tcrvl, "TCRVL");
    strcpy(tcdlt, "TCDLT");
    sprintf(columnnumber, "%d", eventlistfile.cskyx);
    if (fits_read_key(eventlistfile.generic.fptr, TDOUBLE, strcat(tcrpx, columnnumber),
		      &tcrpxx, comment, &status)) break;
    if (fits_read_key(eventlistfile.generic.fptr, TDOUBLE, strcat(tcrvl, columnnumber),
		      &tcrvlx, comment, &status)) break;
    if (fits_read_key(eventlistfile.generic.fptr, TDOUBLE, strcat(tcdlt, columnnumber),
		      &tcdltx, comment, &status)) break;

    // For the sky y-coordinate:
    strcpy(tcrpx, "TCRPX");
    strcpy(tcrvl, "TCRVL");
    strcpy(tcdlt, "TCDLT");
    sprintf(columnnumber, "%d", eventlistfile.cskyy);
    if (fits_read_key(eventlistfile.generic.fptr, TDOUBLE, strcat(tcrpx, columnnumber),
		      &tcrpxy, comment, &status)) break;
    if (fits_read_key(eventlistfile.generic.fptr, TDOUBLE, strcat(tcrvl, columnnumber),
		      &tcrvly, comment, &status)) break;
    if (fits_read_key(eventlistfile.generic.fptr, TDOUBLE, strcat(tcdlt, columnnumber),
		      &tcdlty, comment, &status)) break;
    

    // Get the satellite catalog with the telescope attitude data:
    if (NULL==(attitudecatalog=get_AttitudeCatalog(parameters.attitude_filename,
						   parameters.t0, parameters.timespan, 
						   &status))) break;
						       
    // --- END of Initialization ---



    // --- Beginning of the Sky Imaging Process ---

    // Beginning of actual simulation (after loading required data):
    headas_chat(5, "start sky imaging process ...\n");

    // LOOP over all event in the FITS table
    long attitude_counter=0;   // counter for entries in the AttitudeCatalog

    // SCAN EVENT LIST
    while((EXIT_SUCCESS==status) && (0==EventFileEOF(&eventlistfile.generic))) {
      
      // Read the next event from the FITS file.
      eROSITAEvent event;
      status=eROSITAEventFile_getNextRow(&eventlistfile, &event);
      if(EXIT_SUCCESS!=status) break;

      // Check whether we are finished.
      if (event.time > parameters.t0+parameters.timespan) break;

      // Get the last attitude entry before 'event.time'
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
      // Exact position on the detector:
      struct Point2d detector_position;
      detector_position.x = 
	((double)(event.xi-384/2)+get_random_number()); // in [floating point pixel]
      detector_position.y = 
	((double)(event.yi-384/2)+get_random_number()); // in [floating point pixel]
      double d = sqrt(pow(detector_position.x,2.)+pow(detector_position.y,2.));

      // Determine the off-axis angle corresponding to the detector position.
      double offaxis_angle = d * radperpixel; // [rad]

      // Determine the source position on the sky using the telescope 
      // axis pointing vector and a vector from the point of the intersection 
      // of the optical axis with the sky plane.
      double r = tan(offaxis_angle); // length of this vector

      Vector source_position;
      source_position.x = telescope.nz.x 
	-r*(detector_position.x/d*telescope.nx.x+detector_position.y/d*telescope.ny.x);
      source_position.y = telescope.nz.y 
	-r*(detector_position.x/d*telescope.nx.y+detector_position.y/d*telescope.ny.y);
      source_position.z = telescope.nz.z 
	-r*(detector_position.x/d*telescope.nx.z+detector_position.y/d*telescope.ny.z);
      source_position = normalize_vector(source_position);

      // Determine the equatorial coordinates RA and DEC:
      calculate_ra_dec(source_position, &event.ra, &event.dec);
      event.ra  *= 180./M_PI; // [rad] -> [deg]
      event.dec *= 180./M_PI; // [rad] -> [deg]

      // Determine the pixel coordinates in the sky image:
      event.sky_xi =  (int)((event.ra -tcrvlx)/tcdltx+tcrpxx);
      if ((event.ra -tcrvlx)<0.) {
	event.sky_xi--;
      }
      event.sky_yi =  (int)((event.dec-tcrvly)/tcdlty+tcrpxy);
      if ((event.dec-tcrvly)<0.) {
	event.sky_yi --;
      }

      // Store the data in the Event List FITS file.
      fits_write_col(eventlistfile.generic.fptr, TDOUBLE, eventlistfile.cra,
		     eventlistfile.generic.row, 1, 1, &event.ra, &status);
      fits_write_col(eventlistfile.generic.fptr, TDOUBLE, eventlistfile.cdec, 
		     eventlistfile.generic.row, 1, 1, &event.dec, &status);
      fits_write_col(eventlistfile.generic.fptr, TLONG, eventlistfile.cskyx, 
		     eventlistfile.generic.row, 1, 1, &event.sky_xi, &status);
      fits_write_col(eventlistfile.generic.fptr, TLONG, eventlistfile.cskyy, 
		     eventlistfile.generic.row, 1, 1, &event.sky_yi, &status);

    } // END of scanning-LOOP over the event list.
  } while(0);  // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // release HEADAS random number generator
  HDmtFree();

  // Close the FITS files.
  closeeROSITAEventFile(&eventlistfile);

  // Release memory of AttitudeCatalog
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

  // Get the filename of the input event list (FITS file)
  if ((status = PILGetFname("eventlist_filename", parameters->eventlist_filename))) {
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

  // get the start time of the simulation
  else if ((status = PILGetReal("t0", &parameters->t0))) {
    sprintf(msg, "Error reading the 't0' parameter!\n");
    HD_ERROR_THROW(msg,status);
  }

  // get the timespan for the simulation
  else if ((status = PILGetReal("timespan", &parameters->timespan))) {
    sprintf(msg, "Error reading the 'timespan' parameter!\n");
    HD_ERROR_THROW(msg,status);
  }

  // convert angles from [arc min] to [rad]
  parameters->fov_diameter = parameters->fov_diameter*M_PI/(60.*180.); 

  return(status);
}


