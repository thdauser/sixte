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


struct Parameters {
  char orbit_filename[FILENAME_LENGTH];     // filename of orbit file
  char attitude_filename[FILENAME_LENGTH];  // filename of the attitude file
  char eventlist_filename[FILENAME_LENGTH]; // input: photon list
  char skyimage_filename[FILENAME_LENGTH];  // output: impact list

  double focal_length;
  double fov_diameter;
};


int event_projection_getpar(struct Paramters *parameters);



////////////////////////////////////
// Main procedure.
int event_projection_main() {
  struct Parameters parameters;   // Program parameters

  long sat_nentries; // number of entries in the orbit array ( <= orbit_nrows )
  struct Telescope *sat_catalog=NULL; // catalog with orbit and attitude data 
                                      // over a certain timespan
  struct Eventlist_File eventlistfile;

  fitsfile *impactlist_fptr=NULL;           
  struct Telescope telescope; // Telescope data (like FOV diameter or focal length)

  char msg[MAXMSG];             // error output buffer
  int status=EXIT_SUCCESS;      // error status


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


    // Calculate the minimum cos-value for sources inside the FOV: 
    // (angle(x0,source) <= 1/2 * diameter)
    const double fov_min_align = cos(telescope.fov_diameter/2.); 
    

    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);

    // Open the FITS file with the input event list:
    eventlistfile=open_EventlistFile(parameters.eventlist_filename, &status);
    if ((EXIT_SUCCESS!=status)||(NULL==eventlistfile)) break;

    // Get the satellite catalog with the orbit and (telescope) attitude data:
    if ((status=get_satellite_catalog(&sat_catalog, &sat_nentries, t0, timespan+100., 
				      orbit_filename, attitude_filename))
	!=EXIT_SUCCESS) break;


    // Get the PSF:
    psf = get_psf(psf_filename, &status);
    if (status != EXIT_SUCCESS) break;

    // Create a new FITS file for the output of the impact list:
    remove(impactlist_filename);
    if ((create_impactlist_file(&impactlist_fptr, impactlist_filename, &status))) 
      break;

    
    // --- END of Initialization ---


    // --- Beginning of Imaging Process ---

    // LOOP over all timesteps given the specified timespan from t0 to t0+timespan
    long photonlist_row=0;    // current row in the input list
    long impactlist_row=0;    //      -"-           output list
    long sat_counter=0;       // counter for orbit readout loop


    // Beginning of actual simulation (after loading required data):
    headas_chat(5, "start imaging process ...\n");


    // SCAN PHOTON LIST    
    for(photonlist_row=0; (photonlist_row<photonlist_nrows)&&(status==EXIT_SUCCESS); 
	photonlist_row++) {
      
      // Read an entry from the photon list:
      int anynul = 0;
      struct Photon photon;
      photon.time = 0.;
      photon.energy = 0.;
      photon.ra = 0.;
      photon.dec = 0.;
      fits_read_col(photonlist_fptr, TDOUBLE, 1, photonlist_row+1, 1, 1, 
		    &photon.time, &photon.time, &anynul, &status);
      fits_read_col(photonlist_fptr, TFLOAT, 2, photonlist_row+1, 1, 1, 
		    &photon.energy, &photon.energy, &anynul, &status);
      fits_read_col(photonlist_fptr, TDOUBLE, 3, photonlist_row+1, 1, 1, 
		    &photon.ra, &photon.ra, &anynul, &status);
      fits_read_col(photonlist_fptr, TDOUBLE, 4, photonlist_row+1, 1, 1, 
		    &photon.dec, &photon.dec, &anynul, &status);
      if (status!=EXIT_SUCCESS) break;

      // Rescale from [deg] -> [rad]
      photon.ra  = photon.ra *M_PI/180.;
      photon.dec = photon.dec*M_PI/180.;
      // Determine a unit vector pointing in the direction of the photon.
      photon.direction = unit_vector(photon.ra, photon.dec);


      // Get the last orbit entry before 'photon.time'
      // (in order to interpolate the position and velocity at this time  between 
      // the neighboring calculated orbit positions):
      for( ; sat_counter<sat_nentries; sat_counter++) {
	if(sat_catalog[sat_counter].time>photon.time) {
	  break;
	}
      }
      if(fabs(sat_catalog[--sat_counter].time-photon.time)>600.) { 
	// no entry within 10 minutes !!
	status = EXIT_FAILURE;
	sprintf(msg, "Error: no adequate orbit entry for time %lf!\n", photon.time);
	HD_ERROR_THROW(msg,status);
	break;
      }

      // Check whether the photon is inside the FOV:
      // First determine telescope pointing direction at the actual time.
      telescope.nz = 
	normalize_vector(interpolate_vec(sat_catalog[sat_counter].nz, 
					 sat_catalog[sat_counter].time, 
					 sat_catalog[sat_counter+1].nz, 
					 sat_catalog[sat_counter+1].time, 
					 photon.time));

      // Compare the photon direction to the unit vector specifiing the 
      // direction of the telescope axis:
      if (check_fov(&photon.direction, &telescope.nz, fov_min_align)==0) {
	// Photon is inside the FOV!
	
	// Determine telescope data like direction etc. (attitude).
	// The telescope coordinate system consists of a nx, ny, and nz axis.
	// The nz axis is perpendicular to the detector plane and pointing along
	// the telescope direction. The nx axis is align along the detector 
	// x-direction, which is identical to the detector COLUMN.
	// The ny axis ix pointing along the y-direction of the detector,
	// which is also referred to as ROW.

	// Determine the current nx: perpendicular to telescope axis nz
	// and in the direction of the satellite motion.
	telescope.nx = 
	  normalize_vector(interpolate_vec(sat_catalog[sat_counter].nx, 
					   sat_catalog[sat_counter].time, 
					   sat_catalog[sat_counter+1].nx, 
					   sat_catalog[sat_counter+1].time, 
					   photon.time));
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
	
	// Determine the photon impact position on the detector (in [m]):
	struct Point2d position;  

	// Convolution with PSF:
	// Function returns 0, if the photon does not fall on the detector. 
	// If it hits the detector, the return value is 1.
	if (get_psf_pos(&position, photon, telescope, psf)) {
	  // Check whether the photon hits the detector within the FOV. 
	  // (Due to the effects of the mirrors it might have been scattered over 
	  // the edge of the FOV, although the source is inside the FOV.)
	  if (sqrt(pow(position.x,2.)+pow(position.y,2.)) < 
	      tan(telescope.fov_diameter)*telescope.focal_length) {
	    
	    // Insert the impact position with the photon data into the impact list:
	    fits_insert_rows(impactlist_fptr, impactlist_row++, 1, &status);
	    fits_write_col(impactlist_fptr, TDOUBLE, 1, impactlist_row, 1, 1, 
			   &photon.time, &status);
	    fits_write_col(impactlist_fptr, TFLOAT, 2, impactlist_row, 1, 1, 
			   &photon.energy, &status);
	    fits_write_col(impactlist_fptr, TDOUBLE, 3, impactlist_row, 1, 1, 
			   &position.x, &status);
	    fits_write_col(impactlist_fptr, TDOUBLE, 4, impactlist_row, 1, 1, 
			   &position.y, &status);

	  }
	} // END get_psf_pos(...)
      } // End of FOV check
    } // END of scanning LOOP over the photon list.
  } while(0);  // END of the error handling loop.



  // --- cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // release HEADAS random number generator
  HDmtFree();

  // Close the FITS files.
  if (impactlist_fptr) fits_close_file(impactlist_fptr, &status);
  if (photonlist_fptr) fits_close_file(photonlist_fptr, &status);

  // Release memory of orbit catalog
  if (sat_catalog) free(sat_catalog);

  // Release memory of PSF:
  free_psf(psf);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully!\n\n");

  return(status);
}




////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int event_projection_getpar(struct Paramters *parameters)
{
  char msg[MAXMSG];             // error output buffer
  int status=EXIT_SUCCESS;      // error status

  // Get the filename of the orbit file (FITS file)
  if ((status = PILGetFname("orbit_filename", parameters->orbit_filename))) {
    sprintf(msg, "Error reading the filename of the orbit file!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the filename of the attitude file (FITS file)
  else if ((status = PILGetFname("attitude_filename", parameters->attitude_filename))) {
    sprintf(msg, "Error reading the filename of the attitude file!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the filename of the input photon list (FITS file)
  else if ((status = PILGetFname("eventlist_filename", parameters->eventlist_filename))) {
    sprintf(msg, "Error reading the filename of the event list!\n");
    HD_ERROR_THROW(msg,status);
  }
  
  // Get the filename of the output sky image FITS file
  else if ((status = PILGetFname("skyimage_filename", parameters->skyimage_filename))) {
    sprintf(msg, "Error reading the filename of the sky image output file!\n");
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



