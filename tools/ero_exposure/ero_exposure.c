#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include "sixt.h"
#include "vector.h"
#include "telescope.h"
#include "attitudecatalog.h"
#include "check_fov.h"
#include "vignetting.h"

#define TOOLSUB ero_exposure_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char Attitude[MAXFILENAME];    // filename of the attitude file
  char Vignetting[MAXFILENAME];  // filename of the vignetting file
  char Exposuremap[MAXFILENAME]; // output: exposure map
  
  double t0;
  double timespan;
  /** Step width for the exposure map calculation. */
  double dt; 

  double fov_diameter;

  /** Right ascension range [rad]. */
  double ra1 , ra2;  
  /** Declination range [rad]. */
  double dec1, dec2; 
  /** Number of pixels in right ascension and declination. */
  long ra_bins, dec_bins; 

  /** Projection method (1: AIT, 2: SIN). */
  int projection;

  /** Number of interim maps to be stored. */
  int intermaps;
};


int ero_exposure_getpar(struct Parameters *parameters);


void saveExpoMap(float** const map,
		 const char* const filename,
		 const long naxis1, const long naxis2,
		 struct wcsprm* const wcs,
		 int* const status) 
{
  // 1d image buffer for storing in FITS image.
  float* map1d=NULL;

  // FITS file pointer for exposure map image.
  fitsfile* fptr=NULL;

  // String buffer for FITS header.
  char* headerstr=NULL;

  // Store the exposure map in a FITS file image.
  headas_chat(4, "\nstore exposure map in FITS image '%s' ...\n", 
	      filename);

  do { // Beginning of error handling loop.

    // Convert the exposure map to a 1d-array to store it in the FITS image.
    map1d=(float*)malloc(naxis1*naxis2*sizeof(float));
    CHECK_NULL_BREAK(map1d, *status, 
		    "memory allocation for 1d exposure map buffer failed");
    long x;
    for (x=0; x<naxis1; x++) {
      long y;
      for (y=0; y<naxis2; y++) {
	map1d[x + y*naxis1] = map[x][y];
      }
    }

    // Create a new FITS-file (remove existing one before):
    remove(filename);
    fits_create_file(&fptr, filename, status);
    CHECK_STATUS_BREAK(*status);

    // Create an image in the FITS-file (primary HDU):
    long naxes[2] = { naxis1, naxis2 };
    fits_create_img(fptr, FLOAT_IMG, 2, naxes, status);
    //                               |-> naxis
    CHECK_STATUS_BREAK(*status);

    // Write WCS header keywords.
    int nkeyrec;
    if (0!=wcshdo(0, wcs, &nkeyrec, &headerstr)) {
      SIXT_ERROR("construction of WCS header failed");
      *status=EXIT_FAILURE;
      break;
    }
    char* strptr=headerstr;
    while (strlen(strptr)>0) {
      char strbuffer[81];
      strncpy(strbuffer, strptr, 80);
      strbuffer[80] = '\0';
      fits_write_record(fptr, strbuffer, status);
      CHECK_STATUS_BREAK(*status);
      strptr += 80;
    }
    CHECK_STATUS_BREAK(*status);

    // Write the image to the file.
    long fpixel[2] = {1, 1}; // Lower left corner.
    //                |--|--> FITS coordinates start at (1,1), NOT (0,0).
    // Upper right corner.
    long lpixel[2] = {naxis1, naxis2}; 
    fits_write_subset(fptr, TFLOAT, fpixel, lpixel, map1d, status);
    CHECK_STATUS_BREAK(*status);

  } while(0); // End of error handling loop.

  // Close the exposure map FITS file.
  if(NULL!=fptr) fits_close_file(fptr, status);

  // Release memory.
  if (NULL!=map1d) {
    free(map1d);
  }
  if (NULL!=headerstr) {
    free(headerstr);
  }
}  


int ero_exposure_main() 
{
  // Program parameters.
  struct Parameters par;

  AttitudeCatalog* ac=NULL;
  Vignetting* vignetting=NULL; 
  
  // Array for the calculation of the exposure map.
  float** expoMap=NULL;

  // WCS data structure used for projection.
  struct wcsprm wcs = { .flag=-1 };

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL:
  set_toolname("ero_exposure");
  set_toolversion("0.04");
  

  do {  // Beginning of the ERROR handling loop.

    // --- Initialization ---
    // Read the program parameters using PIL library.
    if ((status=ero_exposure_getpar(&par))) break;

    // Get memory for the exposure map.
    expoMap = (float**)malloc(par.ra_bins*sizeof(float*));
    if (NULL!=expoMap) {
      long x;
      for (x=0; x<par.ra_bins; x++) {
	expoMap[x] = (float*)malloc(par.dec_bins*sizeof(float));
	if (NULL!=expoMap[x]) {
	  // Clear the exposure map.
	  long y;
	  for (y=0; y<par.dec_bins; y++) {
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

    // Set up the WCS data structure.
    if (0!=wcsini(1, 2, &wcs)) {
      SIXT_ERROR("initalization of WCS data structure failed");
      status=EXIT_FAILURE;
      break;
    }
    wcs.naxis =  2;
    wcs.crpix[0] = par.ra_bins/2  + 0.5;
    wcs.crpix[1] = par.dec_bins/2 + 0.5;
    wcs.crval[0] = 0.5*(par.ra1 +par.ra2 )*180./M_PI;
    wcs.crval[1] = 0.5*(par.dec1+par.dec2)*180./M_PI;    
    wcs.cdelt[0] = (par.ra2 -par.ra1 )*180./M_PI/par.ra_bins;
    wcs.cdelt[1] = (par.dec2-par.dec1)*180./M_PI/par.dec_bins;
    strcpy(wcs.cunit[0], "deg");
    strcpy(wcs.cunit[1], "deg");
    if (1==par.projection) {
      strcpy(wcs.ctype[0], "RA---AIT");
      strcpy(wcs.ctype[1], "DEC--AIT");
    } else if (2==par.projection) {
      strcpy(wcs.ctype[0], "RA---SIN");
      strcpy(wcs.ctype[1], "DEC--SIN");
    } else {
      SIXT_ERROR("projection type not supported");
      status=EXIT_FAILURE;
      break;
    }

    // Calculate the minimum cos-value for sources inside the FOV: 
    // (angle(x0,source) <= 1/2 * diameter)
    const double fov_min_align = cos(par.fov_diameter/2.); 
    double field_min_align;
    if ((par.ra2-par.ra1 > M_PI/6.) || 
	(par.dec2-par.dec1 > M_PI/6.)) {
      // Actually -1 should be sufficient, but -2 is even safer.
      field_min_align = -2.; 
    } else {
      field_min_align = 
	cos((sqrt(pow(par.ra2-par.ra1, 2.)+pow(par.dec2-par.dec1, 2.))+
	     par.fov_diameter)/2.);
    }

    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);

    // Get the telescope attitude data.
    ac=loadAttitudeCatalog(par.Attitude, &status);
    CHECK_STATUS_BREAK(status);

    // Load the Vignetting data.
    if (0<strlen(par.Vignetting)) {
      vignetting=newVignetting(par.Vignetting, &status);
      CHECK_STATUS_BREAK(status);
    }

    // --- END of Initialization ---


    // --- Beginning of Exposure Map calculation
    headas_chat(3, "calculate the exposure map ...\n");

    // LOOP over the given time interval from t0 to t0+timespan in steps of dt.
    int intermaps=0;
    double time;
    for (time=par.t0; time<par.t0+par.timespan; time+=par.dt) {
      
      // Print the current time (program status information for the user).
      headas_printf("\rtime: %.1lf s ", time);
      fflush(NULL);

      // Determine the telescope pointing direction at the current time.
      Vector telescope_nz = getTelescopeNz(ac, time, &status);
      CHECK_STATUS_BREAK(status);

      // Calculate the RA and DEC of the pointing direction.
      double telescope_ra, telescope_dec;
      calculate_ra_dec(telescope_nz, &telescope_ra, &telescope_dec);

      // Check if the specified field of the sky might be within the FOV.
      // Otherwise break this run and continue at the beginning of the loop 
      // with the next time step.
      Vector pixpos = 
	unit_vector(0.5*(par.ra1+par.ra2), 0.5*(par.dec1+par.dec2));
      if (check_fov(&pixpos, &telescope_nz, field_min_align)!=0) {
	continue;
      }

      // 2d Loop over the exposure map in order to determine all pixels that
      // are currently within the FOV.
      long x;
      for (x=0; x<par.ra_bins; x++) {
	long y;
	for (y=0; y<par.dec_bins; y++) {
	  
	  double pixcrd[2] = { x+1., y+1. };
	  double imgcrd[2], world[2];
	  double phi, theta;
	  int status2=0;
	  wcsp2s(&wcs, 1, 2, pixcrd, imgcrd, &phi, &theta, world, &status2);
	  if (3==status2) { 
	    // Pixel does not correspond to valid world coordinates.
	    continue;
	  } else if (0!=status2) {
	    SIXT_ERROR("projection failed");
	    status=EXIT_FAILURE;
	    break;
	  }


	  // Determine a unit vector for the calculated RA and Dec.
	  pixpos=unit_vector(world[0]*M_PI/180., world[1]*M_PI/180.);
	  
	  // Check if the current pixel lies within the FOV.
	  if (check_fov(&pixpos, &telescope_nz, fov_min_align)==0) {
	    // Pixel lies inside the FOV!
	
	    // Calculate the off-axis angle ([rad])
	    double delta = acos(scalar_product(&telescope_nz, &pixpos));
	
	    // Add the exposure time step weighted with the vignetting
	    // factor for this particular off-axis angle at 1 keV.
	    double addvalue = par.dt;
	    if (NULL!=vignetting) {
	      addvalue *= get_Vignetting_Factor(vignetting, 1., delta, 0.);
	    }
	    expoMap[x][y] += addvalue;
	  }
	}
	CHECK_STATUS_BREAK(status);  
      }
      CHECK_STATUS_BREAK(status);

      // Check if an interim map should be saved now.
      if (par.intermaps>0) {
	if (time > par.t0+ intermaps*(par.timespan/(par.intermaps+1))) {
	  // Construct the filename.
	  char filename[MAXFILENAME];
	  strncpy(filename, par.Exposuremap, 
		  strlen(par.Exposuremap)-5);
	  filename[strlen(par.Exposuremap)-5]='\0';
	  char buffer[MAXFILENAME];
	  sprintf(buffer, "_%d.fits", intermaps);
	  strcat(filename, buffer);

	  // Save the interim map.
	  saveExpoMap(expoMap, filename, par.ra_bins, par.dec_bins, 
		      &wcs, &status);
	  CHECK_STATUS_BREAK(status);

	  intermaps++;
	}
      }
      // END of saving an interim map.
    } 
    CHECK_STATUS_BREAK(status);
    // END of LOOP over the specified time interval.
    
    // END of generating the exposure map.

    // Store the exposure map in the output file.
    saveExpoMap(expoMap, par.Exposuremap,
		par.ra_bins, par.dec_bins, &wcs, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Release memory.
  freeAttitudeCatalog(&ac);
  destroyVignetting(&vignetting);
  wcsfree(&wcs);

  // Release memory of exposure map.
  if (NULL!=expoMap) {
    long x;
    for (x=0; x<par.ra_bins; x++) {
      if (NULL!=expoMap[x]) {
	free(expoMap[x]);
	expoMap[x]=NULL;
      }
    }
    free(expoMap);
  }

  if (EXIT_SUCCESS==status) headas_chat(2, "finished successfully!\n\n");
  return(status);
}


int ero_exposure_getpar(struct Parameters *par)
{
  int ra_bins, dec_bins;    // Buffer
  int status=EXIT_SUCCESS;  // Error status
  
  // Get the filename of the input attitude file (FITS file)
  if ((status = PILGetFname("Attitude", par->Attitude))) {
    HD_ERROR_THROW("Error reading the filename of the attitude file!\n", status);
  }
  
  // Get the filename of the vignetting data file (FITS file)
  else if ((status = PILGetString("Vignetting", 
				  par->Vignetting))) {
    HD_ERROR_THROW("Error reading the filename of the vignetting file!\n", status);
  }

  // Get the filename of the output exposure map (FITS file)
  else if ((status = PILGetFname("Exposuremap", par->Exposuremap))) {
    HD_ERROR_THROW("Error reading the filename of the exposure map!\n", status);
  }

  // Read the diameter of the FOV (in arcmin)
  else if ((status = PILGetReal("fov_diameter", &par->fov_diameter))) {
    HD_ERROR_THROW("Error reading the diameter of the FOV!\n", status);
  }

  // Get the start time of the exposure map calculation
  else if ((status = PILGetReal("t0", &par->t0))) {
    HD_ERROR_THROW("Error reading the 't0' parameter!\n", status);
  }

  // Get the timespan for the exposure map calculation
  else if ((status = PILGetReal("timespan", &par->timespan))) {
    HD_ERROR_THROW("Error reading the 'timespan' parameter!\n", status);
  }

  // Get the time step for the exposure map calculation
  else if ((status = PILGetReal("dt", &par->dt))) {
    HD_ERROR_THROW("Error reading the 'dt' parameter!\n", status);
  }

  // Get the position of the desired section of the sky 
  // (right ascension and declination range).
  else if ((status = PILGetReal("ra1", &par->ra1))) {
    HD_ERROR_THROW("Error reading the 'ra1' parameter!\n", status);
  }
  else if ((status = PILGetReal("ra2", &par->ra2))) {
    HD_ERROR_THROW("Error reading the 'ra2' parameter!\n", status);
  }
  else if ((status = PILGetReal("dec1", &par->dec1))) {
    HD_ERROR_THROW("Error reading the 'dec1' parameter!\n", status);
  }
  else if ((status = PILGetReal("dec2", &par->dec2))) {
    HD_ERROR_THROW("Error reading the 'dec2' parameter!\n", status);
  }
  // Get the number of x- and y-bins for the exposure map.
  else if ((status = PILGetInt("ra_bins", &ra_bins))) {
    HD_ERROR_THROW("Error reading the number of RA bins!\n", status);
  }
  else if ((status = PILGetInt("dec_bins", &dec_bins))) {
    HD_ERROR_THROW("Error reading the number of DEC bins!\n", status);
  }

  else if ((status = PILGetInt("projection", &par->projection))) {
    HD_ERROR_THROW("Error reading the projection type!\n", status);
  }
  else if ((status = PILGetInt("intermaps", &par->intermaps))) {
    HD_ERROR_THROW("Error reading the number of inter-maps!\n", status);
  }

  // Convert Integer types to Long.
  par->ra_bins  = (long)ra_bins;
  par->dec_bins = (long)dec_bins;

  // Convert angles from [deg] to [rad].
  par->ra1  *= M_PI/180.;
  par->ra2  *= M_PI/180.;
  par->dec1 *= M_PI/180.;
  par->dec2 *= M_PI/180.;

  // Convert angles from [arc min] to [rad].
  par->fov_diameter *= M_PI/(60.*180.); 
  
  return(status);
}

