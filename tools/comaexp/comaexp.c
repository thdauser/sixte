#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include "sixt.h"

#include <wcslib/wcslib.h>

#include "vector.h"
#include "point.h"
#include "telescope.h"
#include "attitudecatalog.h"
#include "check_fov.h"

#define TOOLSUB comaexp_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char attitude_filename[MAXFILENAME];    // filename of the attitude file
  char fovimage_filename[MAXFILENAME];    // filename of the FoV image file
  char exposuremap_filename[MAXFILENAME]; // output: exposure map

  /** Coordinate system: equatorial (0) or galactic (1). */
  int coordinate_system;

  /** Projection method: Plate carrée (0) or Hammer-Aitoff (1). */
  int projection;

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

  // Array for the calculation of the exposure map.
  float** expMap=NULL;       
  struct ImageParameters expMapPar;
  // Array for pre-calculation of the carteesian coordinate vectors
  // of the individual pixels in the exposure map image.
  Vector** pixelpositions=NULL;

  // Image of the FoV.
  float** fovImg=NULL;
  struct ImageParameters fovImgPar;
  // FoV image projection type: 
  // 0: local tangential system (obsolete)
  // 1: Plate carrée (CAR)
  // 2: Gnomonic (TAN)
  int fov_projection;

  // 1-dimensional image buffer for storing in FITS files.
  float*  imagebuffer1d=NULL;     

  long x, y;               // Counters.
  fitsfile* fptr=NULL;     // FITS file pointer for exposure map image.

  int status=EXIT_SUCCESS; // Error status.


  // Register HEATOOL:
  set_toolname("comaexp");
  set_toolversion("0.01");
  

  do { // Beginning of the ERROR handling loop (will at most be run once)

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


    // Read the FoV image file.
    // Open the file.
    fits_open_file(&fptr, parameters.fovimage_filename, READONLY, &status);
    if (EXIT_SUCCESS!=status) break;

    // Determine the width of the image.
    long naxes[2];
    fits_get_img_size(fptr, 2, naxes, &status);
    if (EXIT_SUCCESS!=status) break;
    fovImgPar.ra_bins  = (int)naxes[0];
    fovImgPar.dec_bins = (int)naxes[1];

    // Read the WCS information.
    char comment[MAXMSG];
    char ctype1[MAXMSG], ctype2[MAXMSG];
    fits_read_key(fptr, TSTRING, "CTYPE1", ctype1, 
		  comment, &status);
    fits_read_key(fptr, TSTRING, "CTYPE2", ctype2, 
		  comment, &status);
    fits_read_key(fptr, TDOUBLE, "CDELT1", &fovImgPar.delt1, 
		  comment, &status);
    fits_read_key(fptr, TDOUBLE, "CDELT2", &fovImgPar.delt2, 
		  comment, &status);
    fits_read_key(fptr, TDOUBLE, "CRPIX1", &fovImgPar.rpix1, 
		  comment, &status);
    fits_read_key(fptr, TDOUBLE, "CRPIX2", &fovImgPar.rpix2, 
		  comment, &status);
    fits_read_key(fptr, TDOUBLE, "CRVAL1", &fovImgPar.rval1, 
		  comment, &status);
    fits_read_key(fptr, TDOUBLE, "CRVAL2", &fovImgPar.rval2, 
		  comment, &status);
    if (EXIT_SUCCESS!=status) break;
    // Convert from [deg] to [rad].
    fovImgPar.delt1 *= M_PI/180.;
    fovImgPar.delt2 *= M_PI/180.;
    fovImgPar.rval1 *= M_PI/180.;
    fovImgPar.rval2 *= M_PI/180.;

    // Determine the projection type of the FoV image.
    if ((strlen(ctype1)>0) || (strlen(ctype2)>0)) {
      if ((0==strcmp(&(ctype1[5]), "CAR")) && 
	  (0==strcmp(&(ctype2[5]), "CAR"))) {
	fov_projection = 1;
      } else if ((0==strcmp(&(ctype1[5]), "TAN")) && 
		 (0==strcmp(&(ctype2[5]), "TAN"))) {
	fov_projection = 2;
	// TODO CRVAL2 ???
	fovImgPar.rval2 = fovImgPar.rval1;
      } else {
	status=EXIT_FAILURE;
	HD_ERROR_THROW("Error: FoV image has unknown projection type!\n",
		       status);
	break;
      }
    } else {
      fov_projection = 0;
    }

    // Output of projection type:
    switch(fov_projection) {
    case 0:
      headas_chat(1, "CTYPE of FoV image: none\n");
      break;
    case 1:
      headas_chat(1, "CTYPE of FoV image: CAR (Plate carrée)\n");
      break;
    case 2:
      headas_chat(1, "CTYPE of FoV image: TAN (Gnomonic)\n");
      break;
    default:
      headas_chat(1, "CTYPE of FoV image: unknown\n");
    }

    // Determine the dimensions of the FoV.
    double sin_ra_max = sin(fovImgPar.rval1
			    +(fovImgPar.ra_bins*1.-fovImgPar.rpix1+0.5)*fovImgPar.delt1);
    double sin_ra_min = sin(fovImgPar.rval1-(fovImgPar.rpix1-0.5)*fovImgPar.delt1);
    double sin_dec_max = sin(fovImgPar.rval2
			     +(fovImgPar.dec_bins*1.-fovImgPar.rpix2+0.5)*fovImgPar.delt2);
    double sin_dec_min = sin(fovImgPar.rval2-(fovImgPar.rpix2-0.5)*fovImgPar.delt2);

    headas_chat(5, "FoV dimensions: from %.1lf deg to %.1lf deg (RA direction)\n", 
		asin(sin_ra_min)*180./M_PI, asin(sin_ra_max)*180./M_PI);
    headas_chat(5, "                and from %.1lf deg to %.1lf deg (Dec direction)\n", 
		asin(sin_dec_min)*180./M_PI, asin(sin_dec_max)*180./M_PI);

    // Read the image from the file.
    imagebuffer1d=(float*)malloc(fovImgPar.ra_bins*fovImgPar.dec_bins*
				 sizeof(float));
    if (NULL==imagebuffer1d) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for image buffer failed!\n", status);
      break;
    }      
    int anynul;
    float null_value=0.;
    long fpixel[2] = {1, 1};   // lower left corner
    //                |--|--> FITS coordinates start at (1,1)
    // upper right corner
    long lpixel[2] = {fovImgPar.ra_bins, fovImgPar.dec_bins};  
    long inc[2] = {1, 1};
    fits_read_subset(fptr, TFLOAT, fpixel, lpixel, inc, &null_value, 
		     imagebuffer1d, &anynul, &status);
    if (EXIT_SUCCESS!=status) break;

    // Convert the 1-dimensional image to 2-d.
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
      for (y=0; y<fovImgPar.dec_bins; y++) {
	// Take care of choosing x- and y-axis properly!
	fovImg[x][y]=imagebuffer1d[y*fovImgPar.ra_bins+x];
      }
    }
    if (EXIT_SUCCESS!=status) break;

    // Release memory.
    free(imagebuffer1d);
    imagebuffer1d=NULL;
	
    // Close the file.
    fits_close_file(fptr, &status);
    if (EXIT_SUCCESS!=status) break;
    fptr=NULL;


    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);

    // Get the satellite catalog with the telescope attitude data:
    if (NULL==(ac=loadAttitudeCatalog(parameters.attitude_filename,
				      &status))) break;

    // Pre-calculate the carteesian coordinate vectors of 
    // the positions of the individual pixels in the exposure map.
    for (x=0; x<expMapPar.ra_bins; x++) {
      for (y=0; y<expMapPar.dec_bins; y++) {
	double pixelra  = (x-(expMapPar.rpix1-1.0))*expMapPar.delt1 + expMapPar.rval1;
	double pixeldec = (y-(expMapPar.rpix2-1.0))*expMapPar.delt2 + expMapPar.rval2;
			
	// Check if the requested projection method is Hammer-Aitoff.
	if (1==parameters.projection) {
	  double lon = pixelra * 180./M_PI;
	  double lat = pixeldec* 180./M_PI;
	  double phi   = 0.;
	  double theta = 0.;

	  struct celprm cel;
	  status = celini(&cel);
	  if (EXIT_SUCCESS!=status) break;
	  cel.flag   = 0;
	  cel.offset = 0;
	  cel.phi0   = 0.;
	  cel.theta0 = 0.;
	  cel.ref[0] = 0.;
	  cel.ref[1] = 0.;
	  cel.ref[2] = 0.;
	  cel.ref[3] = 0.;
	  strcpy(cel.prj.code, "AIT");
	  cel.prj.r0     = 0.;
	  //	  cel.prjprm.pv = 
	  cel.prj.phi0   = 0.;
	  cel.prj.theta0 = 0.;

	  // Transform coordinates from projection plane to 
	  // celestial coordinates.
	  int invalid_coordinates=0;
	  status = celx2s(&cel, 1, 1, 1, 1, &lon, &lat,
			  &phi, &theta, &pixelra, &pixeldec, 
			  &invalid_coordinates);

	  if (0==invalid_coordinates) {
	    pixelra = phi   * M_PI/180.;
	    pixeldec= theta * M_PI/180.;
	  } else {
	    pixelpositions[x][y].z = -1000.;
	    status=EXIT_SUCCESS;
	    continue;
	  }
	  if (EXIT_SUCCESS!=status) break;
	}
	// END of check if the requested projection method 
	// is Hammer-Aitoff.

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
	// END of check if requested coordinate system is galactic.

	// Calculate the carteesian coordinate vector.
	pixelpositions[x][y] = unit_vector(pixelra, pixeldec);
      }
      if (EXIT_SUCCESS!=status) break;
      // END of loop over y.
    }    
    if (EXIT_SUCCESS!=status) break;
    // END of loop over x.

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
      Vector nz = getTelescopeNz(ac, time, &status);
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
	  if (pixelpositions[x][y].z < -100.) continue;
	  if (scalar_product(&pixelpositions[x][y], &nz)<0.) continue;

	  // Distinguish between different FoV image projection
	  // types.
	  if (1==fov_projection) {
	    // CAR (Plate carrée).

	    // Check if the pixel is within the telescope FoV.
	    // Projection along y-axis:
	    double sy = scalar_product(&pixelpositions[x][y], &ny);
	    if ((sy > sin_dec_max) || (sy < sin_dec_min)) continue;
	    
	    // Projection along x-axis:
	    double sx = scalar_product(&pixelpositions[x][y], &nx);
	    if ((sx > sin_ra_max) || (sx < sin_ra_min)) continue;

	    double dec = asin(sy);
	    double ra  = asin(sx/cos(dec));
	    
	    int xi = (int)((ra -fovImgPar.rval1)/fovImgPar.delt1+fovImgPar.rpix1+0.5)-1;
	    int yi = (int)((dec-fovImgPar.rval2)/fovImgPar.delt2+fovImgPar.rpix2+0.5)-1;
	    
	    if ((xi<0) || (xi>=fovImgPar.ra_bins )) continue;
	    if ((yi<0) || (yi>=fovImgPar.dec_bins)) continue;

	    expMap[x][y] += parameters.dt * fovImg[xi][yi];

	  } else if (2==fov_projection) {
	    // TAN (Gnomonic).
	    // Use local tangential system with 2 equivalent and
	    // independent angles.

	    // Angle in right ascension direction:
	    double alpha = asin(scalar_product(&pixelpositions[x][y], &nx));
	    // Angle in declination direction:
	    double beta  = asin(scalar_product(&pixelpositions[x][y], &ny));

	    // Image coordinates:
	    int xi = (int)(tan(alpha)/tan(fovImgPar.delt1) + fovImgPar.rpix1 + 0.5) -1;
	    int yi = (int)(tan(beta) /tan(fovImgPar.delt2) + fovImgPar.rpix2 + 0.5) -1;
	    
	    // Check the limits of the FoV.
	    if ((xi >= 0) && (xi < fovImgPar.ra_bins ) &&
		(yi >= 0) && (yi < fovImgPar.dec_bins)) {
	      expMap[x][y] += parameters.dt * fovImg[xi][yi];
	    }

	  } else if (0==fov_projection) {
	    // No particular projection selected for FoV image.
	    // Use local system with 2 equivalent and
	    // independent angles.

	    // Declination direction:
	    double sin_y = scalar_product(&pixelpositions[x][y], &ny);
	    // Right ascension direction:
	    double sin_x = scalar_product(&pixelpositions[x][y], &nx);
	    // Check the limits of the FoV.
	    if ((sin_y < sin_dec_max) && (sin_y > sin_dec_min) &&
		(sin_x < sin_ra_max ) && (sin_x > sin_ra_min )) {
	      double alpha = asin(sin_x);
	      double beta  = asin(sin_y);
	      int xi = (int)((alpha-fovImgPar.rval1)/fovImgPar.delt1+fovImgPar.rpix1+0.5)-1;
	      int yi = (int)((beta -fovImgPar.rval2)/fovImgPar.delt2+fovImgPar.rpix2+0.5)-1;
	      assert(xi>=0);
	      assert(xi<fovImgPar.ra_bins);
	      assert(yi>=0);
	      assert(yi<fovImgPar.dec_bins);
	      expMap[x][y] += parameters.dt * fovImg[xi][yi];
	    }
	  } 
	  // END of different FoV image projection types.
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
    imagebuffer1d = (float*)malloc(parameters.ra_bins*parameters.dec_bins*sizeof(float));
    if (NULL==imagebuffer1d) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for 1d exposure map failed!\n", status);
      break;
    }
    for (x=0; x<parameters.ra_bins; x++) {
      for (y=0; y<parameters.dec_bins; y++) {
	imagebuffer1d[x + y*parameters.ra_bins] = expMap[x][y];
      }
    }

    // Create a new FITS-file (remove existing one before):
    remove(parameters.exposuremap_filename);
    if (fits_create_file(&fptr, parameters.exposuremap_filename, &status)) break;
    // Create an image in the FITS-file (primary HDU):
    naxes[0] = parameters.ra_bins;
    naxes[1] = parameters.dec_bins;
    if (fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status)) break;
    //                                   |-> naxis

    // Store the name of the FoV map in the exposure map FITS file header.
    if (fits_update_key(fptr, TSTRING, "FOVMAP", parameters.fovimage_filename,
			"", &status)) break;   

    // Write WCS keywords to the FITS header of the newly created image.
    double buffer;
    // Use the appropriate coordinate system: either equatorial or 
    // galactic.
    if (0==parameters.coordinate_system) {
      strcpy(ctype1, "RA---");
      strcpy(ctype2, "DEC--");
    } else if (1==parameters.coordinate_system) {
      strcpy(ctype1, "GLON-");
      strcpy(ctype2, "GLAT-");
    }
    if (0==parameters.projection) {
      strcat(ctype1, "CAR");
      strcat(ctype2, "CAR");
    } else if (1==parameters.projection) {
      strcat(ctype1, "AIT");
      strcat(ctype2, "AIT");
    }
    if (fits_update_key(fptr, TSTRING, "CTYPE1", ctype1, "", &status)) break;   
    if (fits_update_key(fptr, TSTRING, "CTYPE2", ctype2, "", &status)) break;   

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
    fpixel[0] = 1; // Lower left corner.
    fpixel[1] = 1; // FITS coordinates start at (1,1), NOT (0,0).
    // Upper right corner.
    lpixel[0] = expMapPar.ra_bins;
    lpixel[1] = expMapPar.dec_bins; 
    fits_write_subset(fptr, TFLOAT, fpixel, lpixel, imagebuffer1d, &status);

  } while(0);  // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Close the exposure map FITS file.
  if(NULL!=fptr) fits_close_file(fptr, &status);

  // Release memory of the attitude catalog.
  freeAttitudeCatalog(&ac);

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
  
  // Image buffer.
  if (NULL!=imagebuffer1d) {
    free(imagebuffer1d);
    imagebuffer1d = NULL;
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
  
  // Get the filename of the attitude file (FITS file)
  if ((status = PILGetFname("attitude_filename", parameters->attitude_filename))) {
    HD_ERROR_THROW("Error reading the filename of the attitude file!\n", status);
  }

  // Get the filename of the FoV image file (FITS file)
  if ((status = PILGetFname("fovimage_filename", parameters->fovimage_filename))) {
    HD_ERROR_THROW("Error reading the filename of the FoV image file!\n", status);
  }
  
  // Get the filename of the output exposure map (FITS file)
  else if ((status = PILGetFname("exposuremap_filename", parameters->exposuremap_filename))) {
    HD_ERROR_THROW("Error reading the filename of the exposure map!\n", status);
  }

  else if ((status = PILGetInt("coordinate_system", &parameters->coordinate_system))) {
    HD_ERROR_THROW("Error reading the type of the coordinate system!\n", status);
  }

  else if ((status = PILGetInt("projection", &parameters->projection))) {
    HD_ERROR_THROW("Error reading the projection method!\n", status);
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



