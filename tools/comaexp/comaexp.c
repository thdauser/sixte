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

  float** expMap=NULL;       // Array for the calculation of the exposure map.
  float*  expMap1d=NULL;     // 1d exposure map for storing in FITS image.
  struct ImageParameters expMapPar;
  // Array for pre-calculation of the carteesian coordinate vectors
  // of the individual pixels in the exposure map image.
  Vector** pixelpositions=NULL;

  // Image of the FoV.
  float** fovImg=NULL;
  struct ImageParameters fovImgPar;

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

    // Determine the WCS parameters of the exposure map.
    expMapPar.ra_bins  = parameters.ra_bins;
    expMapPar.dec_bins = parameters.dec_bins;
    expMapPar.delt1 = (parameters.ra2 -parameters.ra1 )/parameters.ra_bins;
    expMapPar.delt2 = (parameters.dec2-parameters.dec1)/parameters.dec_bins;
    expMapPar.rpix1 = (parameters.ra_bins /2.)+ 0.5;
    expMapPar.rpix2 = (parameters.dec_bins/2.)+ 0.5;
    expMapPar.rval1 = (parameters.ra1 + (expMapPar.ra_bins /2.)*expMapPar.delt1);
    expMapPar.rval2 = (parameters.dec1+ (expMapPar.dec_bins/2.)*expMapPar.delt2);

    // Get memory for the exposure map.
    expMap = (float**)malloc(expMapPar.ra_bins*sizeof(float*));
    if (NULL!=expMap) {
      for (x=0; x<expMapPar.ra_bins; x++) {
	expMap[x] = (float*)malloc(expMapPar.dec_bins*sizeof(float));
	if (NULL!=expMap[x]) {
	  // Clear the exposure map.
	  for (y=0; y<expMapPar.dec_bins; y++) {
	    expMap[x][y] = 0.;
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

    // Get memory for the pixel positions.
    pixelpositions = (Vector**)malloc(expMapPar.ra_bins*sizeof(Vector*));
    if (NULL==pixelpositions) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for exposure map failed!\n", status);
      break;
    }
    for (x=0; x<expMapPar.ra_bins; x++) {
      pixelpositions[x] = (Vector*)malloc(expMapPar.dec_bins*sizeof(Vector));
      if (NULL==pixelpositions[x]) {
	status = EXIT_FAILURE;
	HD_ERROR_THROW("Error: memory allocation for exposure map failed!\n", status);
	break;
      }
    }
    if (EXIT_SUCCESS!=status) break;

    // Set up the dimensions of the FoV image.
    fovImgPar.ra_bins = 512;
    fovImgPar.rpix1   = 256.5;
    fovImgPar.rval1   = 0.;
    fovImgPar.delt1   = 0.1 * M_PI/180.; // [rad]
    double sin_ra_max = sin(fovImgPar.rval1
			    +(fovImgPar.ra_bins*1.-fovImgPar.rpix1+0.5)*fovImgPar.delt1);
    double sin_ra_min = sin(fovImgPar.rval1-(fovImgPar.rpix1-0.5)*fovImgPar.delt1);

    fovImgPar.dec_bins = 512;
    fovImgPar.rpix2    = 256.5;
    fovImgPar.rval2    = -9.48 *M_PI/180.;
    fovImgPar.delt2    =  0.1  *M_PI/180.; // [rad]
    double sin_dec_max = sin(fovImgPar.rval2
			     +(fovImgPar.dec_bins*1.-fovImgPar.rpix2+0.5)*fovImgPar.delt2);
    double sin_dec_min = sin(fovImgPar.rval2-(fovImgPar.rpix2-0.5)*fovImgPar.delt2);

    headas_chat(5, "FoV dimensions: from %.1lf deg to %.1lf deg (RA direction)\n", 
		asin(sin_ra_min)*180./M_PI, asin(sin_ra_max)*180./M_PI);
    headas_chat(5, "            and from %.1lf deg to %.1lf deg (Dec direction)\n", 
		asin(sin_dec_min)*180./M_PI, asin(sin_dec_max)*180./M_PI);

    // Get memory for the FoV image.
    fovImg = (float**)malloc(fovImgPar.ra_bins*sizeof(float*));
    if (NULL==fovImg) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for FoV image failed!\n", status);
      break;
    }
    for (x=0; x<fovImgPar.ra_bins; x++) {
      fovImg[x] = (float*)malloc(fovImgPar.dec_bins*sizeof(float));
      if (NULL==fovImg[x]) {
	status = EXIT_FAILURE;
	HD_ERROR_THROW("Error: memory allocation for FoV image failed!\n", status);
	break;
      }
    }
    if (EXIT_SUCCESS!=status) break;

    // Initialize the FoV image with zero values.
    for (x=0; x<fovImgPar.ra_bins; x++) {
      for (y=0; y<fovImgPar.dec_bins; y++) {
	fovImg[x][y] = 0.;
      }
    }

    // Set up the FoV image for the MIRAX design.
    for (x=0; x<fovImgPar.ra_bins; x++) {
      for (y=0; y<fovImgPar.dec_bins; y++) {
	double ra = (fovImgPar.rval1+(x*1.-fovImgPar.rpix1+1.)*fovImgPar.delt1)*180./M_PI;
	double dec= (fovImgPar.rval2+(y*1.-fovImgPar.rpix2+1.)*fovImgPar.delt2)*180./M_PI;
	if (((ra>0.605)||(ra<-0.605)) && ((dec>-8.88)||(dec<-10.08))){
	  fovImg[x][y] = 1.;
	}
      }
    }
    

    
    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);

    // Get the satellite catalog with the telescope attitude data:
    if (NULL==(ac=get_AttitudeCatalog(parameters.attitude_filename,
				      parameters.t0, parameters.timespan, 
				      &status))) break;

    // Pre-calculate the carteesian coordinate vectors of 
    // the positions of the individual pixels in the exposure map.
    for (x=0; x<expMapPar.ra_bins; x++) {
      for (y=0; y<expMapPar.dec_bins; y++) {
	double pixelra  = (x-(expMapPar.rpix1-1.0))*expMapPar.delt1 + expMapPar.rval1;
	double pixeldec = (y-(expMapPar.rpix2-1.0))*expMapPar.delt2 + expMapPar.rval2;
			
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

	// Calculate the carteesian coordinate vector.
	pixelpositions[x][y] = unit_vector(pixelra, pixeldec);
      }
    }    

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
      for (x=0; x<expMapPar.ra_bins; x++) {
	for (y=0; y<expMapPar.dec_bins; y++) {	  
	  // If the source position is outside the hemisphere
	  // defined by the telescope pointing direction, we
	  // can continue with the next run.
	  // NOTE: This check is necessary in order to avoid
	  // ambiguous values for the subsequent check.
	  if (scalar_product(&pixelpositions[x][y], &nz)<0.) continue;

	  // Check if the pixel is within the telescope FoV.
	  // Declination:
	  double sin_dec = scalar_product(&pixelpositions[x][y], &ny);
	  // Right ascension:
	  double sin_ra  = scalar_product(&pixelpositions[x][y], &nx);
	  // Check the limits of the FoV.
	  if ((sin_dec < sin_dec_max) && (sin_dec > sin_dec_min) &&
	      (sin_ra  < sin_ra_max ) && (sin_ra  > sin_ra_min )) {
	    double ra = asin(sin_ra);
	    double dec= asin(sin_dec);
	    int xi = (int)((ra -fovImgPar.rval1)/fovImgPar.delt1+fovImgPar.rpix1+0.5)-1;
	    int yi = (int)((dec-fovImgPar.rval2)/fovImgPar.delt2+fovImgPar.rpix2+0.5)-1;
	    assert(xi>=0);
	    assert(xi<fovImgPar.ra_bins);
	    assert(yi>=0);
	    assert(yi<fovImgPar.dec_bins);
	    expMap[x][y] += parameters.dt * fovImg[xi][yi];
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
    expMap1d = (float*)malloc(parameters.ra_bins*parameters.dec_bins*sizeof(float));
    if (NULL==expMap1d) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for 1d exposure map failed!\n", status);
      break;
    }
    for (x=0; x<parameters.ra_bins; x++) {
      for (y=0; y<parameters.dec_bins; y++) {
	expMap1d[x + y*parameters.ra_bins] = expMap[x][y];
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
    // Use the appropriate coordinate system: either equatorial or 
    // galactic.
    if (0==parameters.coordinate_system) {
      if (fits_update_key(fptr, TSTRING, "CTYPE1", "RA---CAR", "", &status)) break;   
      if (fits_update_key(fptr, TSTRING, "CTYPE2", "DEC--CAR", "", &status)) break;   
    } else if (1==parameters.coordinate_system) {
      if (fits_update_key(fptr, TSTRING, "CTYPE1", "GLON-CAR", "", &status)) break;   
      if (fits_update_key(fptr, TSTRING, "CTYPE2", "GLAT-CAR", "", &status)) break;         
    }

    if (fits_update_key(fptr, TSTRING, "CUNIT1", "deg", "", &status)) break;   
    buffer = expMapPar.rval1 * 180./M_PI;
    if (fits_update_key(fptr, TDOUBLE, "CRVAL1", &buffer, "", &status)) break;
    buffer = expMapPar.rpix1;
    if (fits_update_key(fptr, TDOUBLE, "CRPIX1", &buffer, "", &status)) break;
    buffer = expMapPar.delt1 * 180./M_PI;
    if (fits_update_key(fptr, TDOUBLE, "CDELT1", &buffer, "", &status)) break;

    if (fits_update_key(fptr, TSTRING, "CUNIT2", "deg", "", &status)) break;   
    buffer = expMapPar.rval2 * 180./M_PI;
    if (fits_update_key(fptr, TDOUBLE, "CRVAL2", &buffer, "", &status)) break;
    buffer = expMapPar.rpix2;
    if (fits_update_key(fptr, TDOUBLE, "CRPIX2", &buffer, "", &status)) break;
    buffer = expMapPar.delt2 * 180./M_PI;
    if (fits_update_key(fptr, TDOUBLE, "CDELT2", &buffer, "", &status)) break;


    // Write the image to the file:
    long fpixel[2] = {1, 1};  // Lower left corner.
    //                |--|--> FITS coordinates start at (1,1), NOT (0,0).
    // Upper right corner.
    long lpixel[2] = {expMapPar.ra_bins, expMapPar.dec_bins}; 
    fits_write_subset(fptr, TFLOAT, fpixel, lpixel, expMap1d, &status);

  } while(0);  // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Close the exposure map FITS file.
  if(NULL!=fptr) fits_close_file(fptr, &status);

  // Release memory of the attitude catalog.
  free_AttitudeCatalog(ac);

  // Release memory of the FoV image.
  if (NULL!=fovImg) {
    for (x=0; x<fovImgPar.ra_bins; x++) {
      if (NULL!=fovImg[x]) {
	free(fovImg[x]);
      }
    }
    free(fovImg);
    fovImg=NULL;
  }

  // Release memory of exposure map.
  if (NULL!=expMap) {
    for (x=0; x<parameters.ra_bins; x++) {
      if (NULL!=expMap[x]) {
	free(expMap[x]);
      }
    }
    free(expMap);
    expMap=NULL;
  }
  if (NULL!=expMap1d) {
    free(expMap1d);
    expMap1d = NULL;
  }

  // Release memory of pixel positions.
  if (NULL!=pixelpositions) {
    for (x=0; x<parameters.ra_bins; x++) {
      if (NULL!=pixelpositions[x]) {
	free(pixelpositions[x]);
      }
    }
    free(pixelpositions);
    pixelpositions=NULL;
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



