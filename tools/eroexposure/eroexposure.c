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
#include "check_fov.h"

#define TOOLSUB eroexposure_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char attitude_filename[FILENAME_LENGTH];    // filename of the attitude file
  char exposuremap_filename[FILENAME_LENGTH]; // output: exposure map
  
  double t0;
  double timespan;

  double fov_diameter;
};


int eroexposure_getpar(struct Parameters *parameters);


////////////////////////////////////
/** Main procedure. */
int eroexposure_main() {
  struct Parameters parameters; // Program parameters.
  
  AttitudeCatalog* attitudecatalog=NULL;
  struct Telescope telescope; // Telescope data (like FOV diameter or focal length).
  const double dt = 1.;       // Step width for the exposure map calculation.
  
  float** exposureMap=NULL;   // Array for the calculation of the exposure map.
  float*  exposureMap1d=NULL; // 1d exposure map for storing in FITS image.
  const long xwidth=288000;   // Dimensions of exposure map.
  const long ywidth=288000;
  int x, y;                   // Counters.
  fitsfile* fptr=NULL;        // FITS file pointer for exposure map image.

  int status=EXIT_SUCCESS;    // Error status.
  char msg[MAXMSG];           // Error output buffer.


  // Register HEATOOL:
  set_toolname("eroexposure");
  set_toolversion("0.01");
  

  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---
    
    // Read the program parameters using PIL library.
    if ((status=eroexposure_getpar(&parameters))) break;

    // Get memory for the exposure map.
    exposureMap = (float**)malloc(xwidth*sizeof(float*));
    if (NULL!=exposureMap) {
      for (x=0; x<xwidth; x++) {
	exposureMap[x] = (float*)malloc(ywidth*sizeof(float));
	if (NULL!=exposureMap[x]) {
	  // Clear the exposure map.
	  for (y=0; y<ywidth; y++) {
	    exposureMap[x][y] = 0.;
	  }
	} else {
	  status = EXIT_FAILURE;
	  HD_ERROR_THROW("Error: memory allocation for exposure map failed!\n", status);
	  break;
	}
      }
    } else {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for exposure map failed!\n", status);
      break;
    }
    if (EXIT_SUCCESS!=status) break;
    
    // Set up the telescope configuration.
    telescope.fov_diameter = parameters.fov_diameter; // [rad]
    // Calculate the minimum cos-value for sources inside the FOV: 
    // (angle(x0,source) <= 1/2 * diameter)
    const double fov_min_align = cos(telescope.fov_diameter/2.); 

    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);

    // Get the satellite catalog with the telescope attitude data:
    if (NULL==(attitudecatalog=get_AttitudeCatalog(parameters.attitude_filename,
						   parameters.t0, parameters.timespan, 
						   &status))) break;
						       
    // --- END of Initialization ---



    // --- Beginning of Exposure Map calculation
    headas_chat(5, "calculation of exposure map ...\n");

    // LOOP over the given time interval from t0 to t0+timespan in steps of dt.
    double time;
    long attitude_counter=0;   // Counter for entries in the AttitudeCatalog.

    for (time=parameters.t0; (time<parameters.t0+parameters.timespan)&&(EXIT_SUCCESS==status);
	 time+=dt) {
      
      // Get the last attitude entry before 'time', in order to interpolate 
      // the attitude at this time between the neighboring calculated values):
      for( ; attitude_counter<attitudecatalog->nentries-1; attitude_counter++) {
	if(attitudecatalog->entry[attitude_counter+1].time > time) {
	  break;
	}
      }
      if(fabs(attitudecatalog->entry[attitude_counter].time-time)>60.) { 
	// No entry within 1 minute !!
	status = EXIT_FAILURE;
	sprintf(msg, "Error: no adequate orbit entry for time %lf!\n", time);
	HD_ERROR_THROW(msg,status);
	break;
      }

      // Determine the telescope pointing direction at the current time.
      telescope.nz = 
	normalize_vector(interpolate_vec(attitudecatalog->entry[attitude_counter].nz, 
					 attitudecatalog->entry[attitude_counter].time, 
					 attitudecatalog->entry[attitude_counter+1].nz, 
					 attitudecatalog->entry[attitude_counter+1].time, 
					 time));
      // Calculate the RA and DEC of the pointing direction.
      double telescope_ra, telescope_dec;
      calculate_ra_dec(telescope.nz, &telescope_ra, &telescope_dec);
      
      // Determine the starting and stopping array indices for the loop over the
      // relevant part of the exposure map:
      double rad_per_pixel = 2*M_PI/1000;
      double x0 = (int)((telescope_ra +M_PI)/rad_per_pixel);
      double y0 = (int)((telescope_dec+M_PI)/rad_per_pixel);
      double x1 = x0 - 100; double x2 = x0 + 100;
      double y1 = y0 - 100; double y2 = y0 + 100;

      // 2d Loop over the exposure map in order to determine all pixels that
      // are currently within the FOV.
      Vector pixel_position;
      for (x=x1; x<=x2; x++) {
	for (y=y1; y<=y2; y++) {
	  pixel_position = unit_vector(telescope_ra  + (x-x0)*rad_per_pixel, 
				       telescope_dec + (y-y0)*rad_per_pixel);
	}
      }

      // Check if the pixel currently lies within the FOV.
      if (check_fov(&pixel_position, &telescope.nz, fov_min_align)==0) {
	// Pixel lies inside the FOV!
	exposureMap[x][y] += dt;
      }
      
    } // END of scanning-LOOP over the specified time interval.
    // END of generating the exposure map.


    // Store the exposure map in a FITS file image.
    headas_chat(5, "store exposure map in FITS image ...\n");

    // Convert the exposure map to a 1d-array to store it in the FITS image.
    exposureMap1d = (float*)malloc(xwidth*ywidth*sizeof(float));
    if (NULL==exposureMap1d) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for 1d exposure map failed!\n", status);
      break;
    }
    for (x=0; x<xwidth; x++) {
      for (y=0; y<ywidth; y++) {
	exposureMap1d[x*xwidth + y] = exposureMap[x][y];
      }
    }

    // Create a new FITS-file:
    if (fits_create_file(&fptr, parameters.exposuremap_filename, &status)) break;
    // Create an image in the FITS-file (primary HDU):
    long naxes[2] = { xwidth, ywidth };
    if (fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status)) break;
    //                                   |-> naxis

    // Write the image to the file:
    long fpixel[2] = {1, 1};  // Lower left corner.
    //                |--|--> FITS coordinates start at (1,1)
    // Upper right corner.
    long lpixel[2] = {xwidth, ywidth}; 
    fits_write_subset(fptr, TFLOAT, fpixel, lpixel, exposureMap1d, &status);

  } while(0);  // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // release HEADAS random number generator
  HDmtFree();

  // Close the exposure map FITS file.
  if(fptr) fits_close_file(fptr, &status);

  // Release memory of AttitudeCatalog
  free_AttitudeCatalog(attitudecatalog);

  // Release memory of exposure map.
  if (NULL!=exposureMap) {
    for (x=0; x<xwidth; x++) {
      if (NULL!=exposureMap[x]) {
	free(exposureMap[x]);
      }
    }
    free(exposureMap);
  }
  if (NULL!=exposureMap1d) {
    free(exposureMap1d);
  }

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully!\n\n");
  return(status);
}




////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int eroexposure_getpar(struct Parameters *parameters)
{
  //  char msg[MAXMSG];         // error output buffer
  int status=EXIT_SUCCESS;  // error status
  
  // Get the filename of the input attitude file (FITS file)
  if ((status = PILGetFname("attitude_filename", parameters->attitude_filename))) {
    HD_ERROR_THROW("Error reading the filename of the attitude file!\n", status);
  }
  
  // Get the filename of the output exposure map (FITS file)
  else if ((status = PILGetFname("exposuremap_filename", parameters->exposuremap_filename))) {
    HD_ERROR_THROW("Error reading the filename of the exposure map!\n", status);
  }

  // Read the diameter of the FOV (in arcmin)
  else if ((status = PILGetReal("fov_diameter", &parameters->fov_diameter))) {
    HD_ERROR_THROW("Error reading the diameter of the FOV!\n", status);
  }

  // Get the start time of the simulation
  else if ((status = PILGetReal("t0", &parameters->t0))) {
    HD_ERROR_THROW("Error reading the 't0' parameter!\n", status);
  }

  // Get the timespan for the simulation
  else if ((status = PILGetReal("timespan", &parameters->timespan))) {
    HD_ERROR_THROW("Error reading the 'timespan' parameter!\n", status);
  }

  // convert angles from [arc min] to [rad]
  parameters->fov_diameter = parameters->fov_diameter*M_PI/(60.*180.); 
  
  return(status);
}



