#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include "sixt.h"
#include "vector.h"
#include "point.h"
#include "telescope.h"
#include "attitudecatalog.h"
#include "check_fov.h"
#include "vignetting.h"

#define TOOLSUB eroexposure_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char attitude_filename[MAXFILENAME];    // filename of the attitude file
  char vignetting_filename[MAXFILENAME];  // filename of the vignetting file
  char exposuremap_filename[MAXFILENAME]; // output: exposure map
  
  double t0;
  double timespan;
  double dt; /**< Step width for the exposure map calculation. */

  double fov_diameter;

  double ra1 , ra2;  /**< Desired right ascension range [rad]. */
  double dec1, dec2; /**< Desired declination range [rad]. */
  long ra_bins, dec_bins; /**< Number of bins in right ascension and declination. */
};


int eroexposure_getpar(struct Parameters *parameters);


struct ImageParameters {
  // WCS parameters.
  // Reference pixels.
  double rpix1, rpix2;
  // [rad] values of reference pixels.
  double rval1, rval2;
  // CDELT1 and CDELT2 WCS-values in [rad].
  double delt1, delt2; 

  long ra_bins, dec_bins;
};
  


static inline int timestep(float** expoMap, 
			   struct ImageParameters* params,
			   double time, double dt,
			   AttitudeCatalog* ac, 
			   Vignetting* vignetting,
			   double field_align, double fov_align)
{
  double telescope_ra, telescope_dec;
  Vector telescope_nz, pixel_position;

  int status = EXIT_SUCCESS;

  // Determine the telescope pointing direction at the current time.
  telescope_nz = getTelescopeNz(ac, time, &status);
  if (EXIT_SUCCESS!=status) return(status);

  // Calculate the RA and DEC of the pointing direction.
  calculate_ra_dec(telescope_nz, &telescope_ra, &telescope_dec);

  // Check if the specified field of the sky might be within the FOV.
  // Otherwise break this run and continue at the beginning of the loop 
  // with the next time step.
  pixel_position = unit_vector((params->ra_bins /2-(params->rpix1-0.5))*params->delt1 + 
			       params->rval1,
			       (params->dec_bins/2-(params->rpix2-0.5))*params->delt2 + 
			       params->rval2);
  if (check_fov(&pixel_position, &telescope_nz, field_align)!=0) return(status);

  // 2d Loop over the exposure map in order to determine all pixels that
  // are currently within the FOV.
  long x, y;                  // Counters.
  long x1, x2, y1, y2;        
  x1 = 0; x2 = params->ra_bins -1;
  y1 = 0; y2 = params->dec_bins-1;
  // Buffer for off-axis angle:
  double theta;
  for (x=x1; x<=x2; x++) {
    for (y=y1; y<=y2; y++) {
      pixel_position = unit_vector((x-(params->rpix1-1.0))*params->delt1 + params->rval1,
				   (y-(params->rpix2-1.0))*params->delt2 + params->rval2);
	    
      // Check if the pixel lies CLOSE to the FOV.
      // If not make a bigger jump.
	  
      // Check if the current pixel lies within the FOV.
      if (check_fov(&pixel_position, &telescope_nz, fov_align)==0) {
	// Pixel lies inside the FOV!
	
	// Calculate the off-axis angle ([rad])
	theta = acos(scalar_product(&telescope_nz, &pixel_position));
	
	// Add the exposure time step weighted with the vignetting
	// factor for this particular off-axis angle at 1 keV.
	// The azimuthal angle is neglected (TODO).
	expoMap[x][y] += 
	  dt* get_Vignetting_Factor(vignetting, 1., theta, 0.);
      }
    }
  }

  return(status);
}


////////////////////////////////////
/** Main procedure. */
int eroexposure_main() {
  struct Parameters parameters; // Program parameters.
  
  AttitudeCatalog* ac=NULL;
  // Mirror vignetting data.
  Vignetting* vignetting=NULL; 
  
  float** expoMap=NULL;       // Array for the calculation of the exposure map.
  float*  expoMap1d=NULL;     // 1d exposure map for storing in FITS image.
  struct ImageParameters imgParams;
  long x, y;                  // Counters.
  fitsfile* fptr=NULL;        // FITS file pointer for exposure map image.

  int status=EXIT_SUCCESS;    // Error status.


  // Register HEATOOL:
  set_toolname("eroexposure");
  set_toolversion("0.01");
  

  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---
    // Read the program parameters using PIL library.
    if ((status=eroexposure_getpar(&parameters))) break;

    imgParams.ra_bins  = parameters.ra_bins;
    imgParams.dec_bins = parameters.dec_bins;

    // Get memory for the exposure map.
    expoMap = (float**)malloc(imgParams.ra_bins*sizeof(float*));
    if (NULL!=expoMap) {
      for (x=0; x<imgParams.ra_bins; x++) {
	expoMap[x] = (float*)malloc(imgParams.dec_bins*sizeof(float));
	if (NULL!=expoMap[x]) {
	  // Clear the exposure map.
	  for (y=0; y<imgParams.dec_bins; y++) {
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
    imgParams.delt1 = (parameters.ra2 -parameters.ra1 )/parameters.ra_bins;
    imgParams.delt2 = (parameters.dec2-parameters.dec1)/parameters.dec_bins;
    imgParams.rpix1 = (parameters.ra_bins /2)+ 0.5;
    imgParams.rpix2 = (parameters.dec_bins/2)+ 0.5;
    imgParams.rval1 = (parameters.ra1 + (imgParams.ra_bins /2)*imgParams.delt1);
    imgParams.rval2 = (parameters.dec1+ (imgParams.dec_bins/2)*imgParams.delt2);
    
    // Set up the telescope configuration.
    float fov_diameter = parameters.fov_diameter; // [rad]
    // Calculate the minimum cos-value for sources inside the FOV: 
    // (angle(x0,source) <= 1/2 * diameter)
    const double fov_min_align = cos(fov_diameter/2.); 
    double field_min_align;
    if ((parameters.ra2-parameters.ra1 > M_PI/6.) || 
	(parameters.dec2-parameters.dec1 > M_PI/6.)) {
      field_min_align = -2.; // Actually -1 should be sufficient, but -2 is even safer.
    } else {
      field_min_align = cos((sqrt(pow(parameters.ra2-parameters.ra1, 2.) +
				  pow(parameters.dec2-parameters.dec1, 2.)) +
			     fov_diameter)/2.);
    }

    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);

    // Get the satellite catalog with the telescope attitude data:
    if (NULL==(ac=loadAttitudeCatalog(parameters.attitude_filename,
				      parameters.t0, parameters.timespan, 
				      &status))) break;

    // Get the Vignetting data:
    vignetting = newVignetting(parameters.vignetting_filename, &status);
    if (status != EXIT_SUCCESS) break;

    // --- END of Initialization ---


    // --- Beginning of Exposure Map calculation
    headas_chat(5, "calculate the exposure map ...\n");

    //#pragma omp parallel
    //    {
      // LOOP over the given time interval from t0 to t0+timespan in steps of dt.
      double time;
      //#pragma omp for private(status)
      for (time=parameters.t0; time<parameters.t0+parameters.timespan;
	   time+=parameters.dt) {
      
	// Print the current time (program status information for the user).
	//headas_printf("\rtime: %.1lf s ", time);
	//fflush(NULL);

	status=timestep(expoMap, &imgParams, time, parameters.dt,
			ac, vignetting, field_min_align, fov_min_align);
	if (status != EXIT_SUCCESS) break;
      } 
      if (status != EXIT_SUCCESS) break;
      // END of LOOP over the specified time interval.
      //    }
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
    if (fits_update_key(fptr, TSTRING, "CTYPE1", "RA---CAR", "", &status)) break;   
    if (fits_update_key(fptr, TSTRING, "CUNIT1", "deg", "", &status)) break;   
    buffer = imgParams.rval1 * 180./M_PI;
    if (fits_update_key(fptr, TDOUBLE, "CRVAL1", &buffer, "", &status)) break;
    buffer = imgParams.rpix1;
    if (fits_update_key(fptr, TDOUBLE, "CRPIX1", &buffer, "", &status)) break;
    buffer = imgParams.delt1 * 180./M_PI;
    if (fits_update_key(fptr, TDOUBLE, "CDELT1", &buffer, "", &status)) break;

    if (fits_update_key(fptr, TSTRING, "CTYPE2", "DEC--CAR", "", &status)) break;   
    if (fits_update_key(fptr, TSTRING, "CUNIT2", "deg", "", &status)) break;   
    buffer = imgParams.rval2 * 180./M_PI;
    if (fits_update_key(fptr, TDOUBLE, "CRVAL2", &buffer, "", &status)) break;
    buffer = imgParams.rpix2;
    if (fits_update_key(fptr, TDOUBLE, "CRPIX2", &buffer, "", &status)) break;
    buffer = imgParams.delt2 * 180./M_PI;
    if (fits_update_key(fptr, TDOUBLE, "CDELT2", &buffer, "", &status)) break;


    // Write the image to the file:
    long fpixel[2] = {1, 1};  // Lower left corner.
    //                |--|--> FITS coordinates start at (1,1), NOT (0,0).
    // Upper right corner.
    long lpixel[2] = {imgParams.ra_bins, imgParams.dec_bins}; 
    fits_write_subset(fptr, TFLOAT, fpixel, lpixel, expoMap1d, &status);

  } while(0);  // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Close the exposure map FITS file.
  if(NULL!=fptr) fits_close_file(fptr, &status);

  // Release memory.
  freeAttitudeCatalog(&ac);
  destroyVignetting(&vignetting);

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
  
  // Get the filename of the vignetting data file (FITS file)
  else if ((status = PILGetFname("vignetting_filename", parameters->vignetting_filename))) {
    HD_ERROR_THROW("Error reading the filename of the vignetting file!\n", status);
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



