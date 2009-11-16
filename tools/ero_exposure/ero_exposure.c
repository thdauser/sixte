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
  double dt; /**< Step width for the exposure map calculation. */

  double fov_diameter;

  double ra1 , ra2;  /**< Desired right ascension range [rad]. */
  double dec1, dec2; /**< Desired declination range [rad]. */
  long ra_bins, dec_bins; /**< Number of bins in right ascension and declination. */
};


int eroexposure_getpar(struct Parameters *parameters);


////////////////////////////////////
/** Main procedure. */
int eroexposure_main() {
  struct Parameters parameters; // Program parameters.
  
  AttitudeCatalog* attitudecatalog=NULL;
  struct Telescope telescope; // Telescope data (like FOV diameter or focal length).
  
  float** expoMap=NULL;       // Array for the calculation of the exposure map.
  float*  expoMap1d=NULL;     // 1d exposure map for storing in FITS image.
  long x, y;                  // Counters.
  fitsfile* fptr=NULL;        // FITS file pointer for exposure map image.

  // WCS parameters.
  // Reference pixels.
  double rpix1, rpix2;
  // [rad] values of reference pixels.
  double rval1, rval2;
  // CDELT1 and CDELT2 WCS-values in [rad].
  double delt1, delt2; 

  double telescope_ra, telescope_dec;
  Vector pixel_position;

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
    expoMap = (float**)malloc(parameters.ra_bins*sizeof(float*));
    if (NULL!=expoMap) {
      for (x=0; x<parameters.ra_bins; x++) {
	expoMap[x] = (float*)malloc(parameters.dec_bins*sizeof(float));
	if (NULL!=expoMap[x]) {
	  // Clear the exposure map.
	  for (y=0; y<parameters.dec_bins; y++) {
	    expoMap[x][y] = 0.;
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

    // Determine the WCS parameters.
    delt1 = (parameters.ra2 -parameters.ra1 )/parameters.ra_bins;
    delt2 = (parameters.dec2-parameters.dec1)/parameters.dec_bins;
    rpix1 = (parameters.ra_bins /2)+ 0.5;
    rpix2 = (parameters.dec_bins/2)+ 0.5;
    rval1 = (parameters.ra1 + (parameters.ra_bins /2)*delt1);
    rval2 = (parameters.dec1+ (parameters.dec_bins/2)*delt2);
    
    // Set up the telescope configuration.
    telescope.fov_diameter = parameters.fov_diameter; // [rad]
    // Calculate the minimum cos-value for sources inside the FOV: 
    // (angle(x0,source) <= 1/2 * diameter)
    const double fov_min_align = cos(telescope.fov_diameter/2.); 
    double field_min_align;
    if ((parameters.ra2-parameters.ra1 > M_PI/6.) || 
	(parameters.dec2-parameters.dec1 > M_PI/6.)) {
      field_min_align = -2.; // Actually -1 should be sufficient, but -2 is even safer.
    } else {
      field_min_align = cos((sqrt(pow(parameters.ra2-parameters.ra1, 2.) +
				  pow(parameters.dec2-parameters.dec1, 2.)) +
			     telescope.fov_diameter)/2.);
    }

    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);

    // Get the satellite catalog with the telescope attitude data:
    if (NULL==(attitudecatalog=get_AttitudeCatalog(parameters.attitude_filename,
						   parameters.t0, parameters.timespan, 
						   &status))) break;
    // --- END of Initialization ---



    // --- Beginning of Exposure Map calculation
    headas_chat(5, "calculate the exposure map ...\n");

    // LOOP over the given time interval from t0 to t0+timespan in steps of dt.
    double time;
    long attitude_counter=0; // Counter for entries in the AttitudeCatalog.

    for (time=parameters.t0; (time<parameters.t0+parameters.timespan)&&(EXIT_SUCCESS==status);
	 time+=parameters.dt) {
      
      // Print the current time (program status information for the user).
      headas_printf("\rtime: %.1lf s ", time);
      fflush(NULL);

      // Get the last attitude entry before 'time', in order to interpolate 
      // the attitude at this time between the neighboring calculated values):
      for( ; attitude_counter<attitudecatalog->nentries-1; attitude_counter++) {
	if(attitudecatalog->entry[attitude_counter+1].time > time) {
	  break;
	}
      }
      if(fabs(attitudecatalog->entry[attitude_counter].time-time)>600.) { 
	// No entry within 1 minute !!
	status = EXIT_FAILURE;
	sprintf(msg, "Error: no adequate orbit entry for time %lf!\n", time);
	HD_ERROR_THROW(msg, status);
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
      calculate_ra_dec(telescope.nz, &telescope_ra, &telescope_dec);

      // Check if the specified field of the sky might be within the FOV.
      // Otherwise break this run and continue at the beginning of the loop 
      // with the next time step.
      pixel_position = unit_vector((parameters.ra_bins/2-(rpix1-1.0))*delt1 + rval1,
				   (parameters.dec_bins/2-(rpix2-1.0))*delt2 + rval2);
      if (check_fov(&pixel_position, &telescope.nz, field_min_align)!=0) continue;

      // 2d Loop over the exposure map in order to determine all pixels that
      // are currently within the FOV.
      double x1 = 0; double x2 = parameters.ra_bins-1;
      double y1 = 0; double y2 = parameters.dec_bins-1;
      for (x=x1; x<=x2; x++) {
	for (y=y1; y<=y2; y++) {
	  pixel_position = unit_vector((x-(rpix1-1.0))*delt1 + rval1,
				       (y-(rpix2-1.0))*delt2 + rval2);
				       
	  // Check if the pixel lies CLOSE to the FOV.
	  // If not make a bigger jump.
	  
	  // Check if the current pixel lies within the FOV.
	  if (check_fov(&pixel_position, &telescope.nz, fov_min_align)==0) {
	    // Pixel lies inside the FOV!
	    long xi=x, yi=y;
	    /*
	    while (xi<                   0) xi+=parameters.ra_bins;
	    while (xi>= parameters.ra_bins) xi-=parameters.ra_bins;
	    while (yi<                   0) yi+=parameters.dec_bins;
	    while (yi>=parameters.dec_bins) yi-=parameters.dec_bins;
	    */
	    expoMap[xi][yi] += parameters.dt;
	  }
	}
      }
    } // END of scanning-LOOP over the specified time interval.
    // END of generating the exposure map.


    // Store the exposure map in a FITS file image.
    headas_chat(5, "\nstore exposure map in FITS image '%s' ...\n", 
		parameters.exposuremap_filename);

    // Convert the exposure map to a 1d-array to store it in the FITS image.
    expoMap1d = (float*)malloc(parameters.ra_bins*parameters.dec_bins*sizeof(float));
    if (NULL==expoMap1d) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for 1d exposure map failed!\n", status);
      break;
    }
    for (x=0; x<parameters.ra_bins; x++) {
      for (y=0; y<parameters.dec_bins; y++) {
	expoMap1d[x + y*parameters.ra_bins] = expoMap[x][y];
      }
    }

    // Create a new FITS-file (remove existing one before):
    remove(parameters.exposuremap_filename);
    if (fits_create_file(&fptr, parameters.exposuremap_filename, &status)) break;
    // Create an image in the FITS-file (primary HDU):
    long naxes[2] = { parameters.ra_bins, parameters.dec_bins };
    if (fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status)) break;
    //                                   |-> naxis

    // Write WCS keywords to the FITS header of the newly created image.
    double buffer;
    if (fits_update_key(fptr, TSTRING, "CTYPE1", "RA", "", &status)) break;   
    if (fits_update_key(fptr, TSTRING, "CUNIT1", "deg", "", &status)) break;   
    buffer = rval1 * 180./M_PI;
    if (fits_update_key(fptr, TDOUBLE, "CRVAL1", &buffer, "", &status)) break;
    buffer = rpix1;
    if (fits_update_key(fptr, TDOUBLE, "CRPIX1", &buffer, "", &status)) break;
    buffer = delt1 * 180./M_PI;
    if (fits_update_key(fptr, TDOUBLE, "CDELT1", &buffer, "", &status)) break;

    if (fits_update_key(fptr, TSTRING, "CTYPE2", "DEC", "", &status)) break;   
    if (fits_update_key(fptr, TSTRING, "CUNIT2", "deg", "", &status)) break;   
    buffer = rval2 * 180./M_PI;
    if (fits_update_key(fptr, TDOUBLE, "CRVAL2", &buffer, "", &status)) break;
    buffer = rpix2;
    if (fits_update_key(fptr, TDOUBLE, "CRPIX2", &buffer, "", &status)) break;
    buffer = delt2  * 180./M_PI;
    if (fits_update_key(fptr, TDOUBLE, "CDELT2", &buffer, "", &status)) break;


    // Write the image to the file:
    long fpixel[2] = {1, 1};  // Lower left corner.
    //                |--|--> FITS coordinates start at (1,1), NOT (0,0).
    // Upper right corner.
    long lpixel[2] = {parameters.ra_bins, parameters.dec_bins}; 
    fits_write_subset(fptr, TFLOAT, fpixel, lpixel, expoMap1d, &status);

  } while(0);  // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Close the exposure map FITS file.
  if(NULL!=fptr) fits_close_file(fptr, &status);

  // Release memory of AttitudeCatalog.
  free_AttitudeCatalog(attitudecatalog);

  // Release memory of exposure map.
  if (NULL!=expoMap) {
    for (x=0; x<parameters.ra_bins; x++) {
      if (NULL!=expoMap[x]) {
	free(expoMap[x]);
	expoMap[x]=NULL;
      }
    }
    free(expoMap);
  }
  if (NULL!=expoMap1d) {
    free(expoMap1d);
  }

  if (EXIT_SUCCESS==status) headas_chat(5, "finished successfully!\n\n");
  return(status);
}




////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int eroexposure_getpar(struct Parameters *parameters)
{
  int ra_bins, dec_bins;    // Buffer
  int status=EXIT_SUCCESS;  // Error status
  
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

  // Get the start time of the exposure map calculation
  else if ((status = PILGetReal("t0", &parameters->t0))) {
    HD_ERROR_THROW("Error reading the 't0' parameter!\n", status);
  }

  // Get the timespan for the exposure map calculation
  else if ((status = PILGetReal("timespan", &parameters->timespan))) {
    HD_ERROR_THROW("Error reading the 'timespan' parameter!\n", status);
  }

  // Get the time step for the exposure map calculation
  else if ((status = PILGetReal("dt", &parameters->dt))) {
    HD_ERROR_THROW("Error reading the 'dt' parameter!\n", status);
  }

  // Get the position of the desired section of the sky 
  // (right ascension and declination range).
  else if ((status = PILGetReal("ra1", &parameters->ra1))) {
    HD_ERROR_THROW("Error reading the 'ra1' parameter!\n", status);
  }
  else if ((status = PILGetReal("ra2", &parameters->ra2))) {
    HD_ERROR_THROW("Error reading the 'ra2' parameter!\n", status);
  }
  else if ((status = PILGetReal("dec1", &parameters->dec1))) {
    HD_ERROR_THROW("Error reading the 'dec1' parameter!\n", status);
  }
  else if ((status = PILGetReal("dec2", &parameters->dec2))) {
    HD_ERROR_THROW("Error reading the 'dec2' parameter!\n", status);
  }
  // Get the number of bins for the exposure map.
  else if ((status = PILGetInt("ra_bins", &ra_bins))) {
    HD_ERROR_THROW("Error reading the number of RA bins!\n", status);
  }
  else if ((status = PILGetInt("dec_bins", &dec_bins))) {
    HD_ERROR_THROW("Error reading the number of DEC bins!\n", status);
  }

  // Convert Integer types to Long.
  parameters->ra_bins  = (long)ra_bins;
  parameters->dec_bins = (long)dec_bins;

  // Convert angles from [deg] to [rad].
  parameters->ra1  *= M_PI/180.;
  parameters->ra2  *= M_PI/180.;
  parameters->dec1 *= M_PI/180.;
  parameters->dec2 *= M_PI/180.;

  // Convert angles from [arc min] to [rad].
  parameters->fov_diameter *= M_PI/(60.*180.); 
  
  return(status);
}



