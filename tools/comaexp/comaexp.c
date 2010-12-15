#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include "sixt.h"
#include "sixt_random.h"
#include "vector.h"
#include "point.h"
#include "telescope.h"
#include "attitudecatalog.h"
#include "check_fov.h"

#define TOOLSUB comaexp_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char attitude_filename[FILENAME_LENGTH];    // filename of the attitude file
  char exposuremap_filename[FILENAME_LENGTH]; // output: exposure map

  /** Coordinate system: equatorial (0) or galactic (1). */
  int coordinate_system;

  double t0;
  double timespan;
  double dt; /**< Step width for the exposure map calculation. */

  double ra1 , ra2;  /**< Desired right ascension range [rad]. */
  double dec1, dec2; /**< Desired declination range [rad]. */
  long ra_bins, dec_bins; /**< Number of bins in right ascension and declination. */
};


int comaexp_getpar(struct Parameters *parameters);


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
  


////////////////////////////////////
/** Main procedure. */
int comaexp_main() 
{
  struct Parameters parameters; // Program parameters.
  
  AttitudeCatalog* ac=NULL;
  
  float** expoMap=NULL;       // Array for the calculation of the exposure map.
  float*  expoMap1d=NULL;     // 1d exposure map for storing in FITS image.
  struct ImageParameters imgParams;
  long x, y;                  // Counters.
  fitsfile* fptr=NULL;        // FITS file pointer for exposure map image.

  int status=EXIT_SUCCESS;    // Error status.


  // Register HEATOOL:
  set_toolname("comaexp");
  set_toolversion("0.01");
  

  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---
    // Read the program parameters using PIL library.
    if ((status=comaexp_getpar(&parameters))) break;

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
    
    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);

    // Get the satellite catalog with the telescope attitude data:
    if (NULL==(ac=get_AttitudeCatalog(parameters.attitude_filename,
				      parameters.t0, parameters.timespan, 
				      &status))) break;

    // --- END of Initialization ---


    // --- Beginning of Exposure Map calculation
    headas_chat(5, "calculate the exposure map ...\n");

    // LOOP over the given time interval from t0 to t0+timespan in steps of dt.
    double time;
    for (time=parameters.t0; time<parameters.t0+parameters.timespan;
	 time+=parameters.dt) {
      
      // Print the current time (program status information for the user).
      headas_printf("\rtime: %.1lf s ", time);
      fflush(NULL);

      // Determine the telescope pointing direction.
      Vector nz = getTelescopePointing(ac, time, &status);
      if (status != EXIT_SUCCESS) break;

      // Determine the two other axes of the telescope reference
      // system (carteesian).
      // Reference vectors:
      Vector north = { .x = 0., .y = 0., .z = 1. };
      Vector n2 = vector_product(normalize_vector(vector_product(nz, north)),nz);
      Vector n1 = vector_product(nz, n2);
      // Consider the roll-angle:
      double roll_angle = getRollAngle(ac, time, &status); // [rad]
      if (status != EXIT_SUCCESS) break;
      Vector nx = {
	.x = n1.x * cos(roll_angle) + n2.x * sin(roll_angle),
	.y = n1.y * cos(roll_angle) + n2.y * sin(roll_angle),
	.z = n1.z * cos(roll_angle) + n2.z * sin(roll_angle)
      };
      Vector ny = {
	.x = - n1.x * sin(roll_angle) + n2.x * cos(roll_angle),
	.y = - n1.y * sin(roll_angle) + n2.y * cos(roll_angle),
	.z = - n1.z * sin(roll_angle) + n2.z * cos(roll_angle)
      };

      // Loop over all pixel in the exposure map.
      for (x=0; x<imgParams.ra_bins; x++) {
	for (y=0; y<imgParams.dec_bins; y++) {	  
	  // Determine the pointing vector to the pixel.
	  double pixelra  = (x-(imgParams.rpix1-1.0))*imgParams.delt1 + imgParams.rval1;
	  double pixeldec = (y-(imgParams.rpix2-1.0))*imgParams.delt2 + imgParams.rval2;
			
	  // If the coordinate system of the exposure map is galactic 
	  // coordinates, convert the pixel position vector from 
	  // galactic to equatorial coordinates.
	  if (1==parameters.coordinate_system) {
	    double lon = pixelra;
	    double lat = pixeldec;
	    const double l_ncp = 2.145566759798267518;
	    const double cos_d_ngp = 0.8899880874849542;
	    const double sin_d_ngp = 0.4559837761750669;
	    pixelra  = 
	      atan2(cos(lat)*sin(l_ncp - lon), 
		    cos_d_ngp*sin(lat)-sin_d_ngp*cos(lat)*cos(l_ncp - lon)) +
	      +3.3660332687500039;
	    while (pixelra>2*M_PI) {
	      pixelra -= 2*M_PI;
	    }
	    while (pixelra<0.) {
	      pixelra += 2*M_PI;
	    }
	    pixeldec = asin(sin_d_ngp*sin(lat) + cos_d_ngp*cos(lat)*cos(l_ncp - lon));
	  }

	  // Calculate a carteesian coordinate vector.
	  Vector pixelpos = unit_vector(pixelra, pixeldec);

	  // If the source position is outside the hemisphere
	  // defined by the telescope pointing direction, we
	  // can continue with the next run.
	  // NOTE: This check is necessary in order to avoid
	  // ambiguous values for the subsequent check.
	  if (scalar_product(&pixelpos, &nz)<0.) continue;

	  // Define the FoV.
	  const double dec_max =  15. *M_PI/180.; // [rad]
	  const double dec_min = -35. *M_PI/180.; // [rad]
	  const double ra_max  =  25. *M_PI/180.; // [rad]
	  const double ra_min  = -25. *M_PI/180.; // [rad]
	  
	  // Check if the pixel is within the telescope FoV.
	  // Declination:
	  double dec = scalar_product(&pixelpos, &ny);
	  // Right ascension:
	  double ra  = scalar_product(&pixelpos, &nx);
	  // Check the limits of the FoV.
	  if ((dec < sin(dec_max)) && (dec > sin(dec_min)) &&
	      (ra  < sin(ra_max) ) && (ra  > sin(ra_min) )) {
	    expoMap[x][y] += parameters.dt;
	  }
	}
      }
      if (status != EXIT_SUCCESS) break;
      // END of loop over all pixels in the exposure map.
    } 
    if (status != EXIT_SUCCESS) break;
    // END of LOOP over the specified time interval.
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
  free_AttitudeCatalog(ac);

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
int comaexp_getpar(struct Parameters *parameters)
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

  else if ((status = PILGetInt("coordinate_system", &parameters->coordinate_system))) {
    HD_ERROR_THROW("Error reading the type of the coordinate system!\n", status);
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

  if ((parameters->coordinate_system<0)||(parameters->coordinate_system>1)) {
    status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: invalid coordinate system!\n", status);
  }

  return(status);
}



