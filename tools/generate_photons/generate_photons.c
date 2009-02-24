#include "generate_photons.h"



////////////////////////////
int generate_photons_getpar(
			    char orbit_filename[],
			    char attitude_filename[],
			    char photonlist_filename[],
			    double *t0,
			    double *timespan,
			    double *bandwidth
			    )
{
  char msg[MAXMSG];           // error message buffer
  int status = EXIT_SUCCESS;  // error status flag

  // Get the filename of the Orbit file (FITS file):
  if ((status = PILGetFname("orbit_filename", orbit_filename))) {
    sprintf(msg, "Error reading the filename of the orbit file!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the filename of the Attitude file (FITS file):
  else if ((status = PILGetFname("attitude_filename", attitude_filename))) {
    sprintf(msg, "Error reading the filename of the attitude file!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the filename of the Photon-List file (FITS file):
  else if ((status = PILGetFname("photonlist_filename", photonlist_filename))) {
    sprintf(msg, "Error reading the filename of the output file for "
	    "the photon list!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the start time of the photon generation
  else if ((status = PILGetReal("t0", t0))) {
    sprintf(msg, "Error reading the 't0' parameter!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the timespan for the photon generation
  else if ((status = PILGetReal("timespan", timespan))) {
    sprintf(msg, "Error reading the 'timespan' parameter!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the (half) width of the preselection band [arcmin]
  else if ((status = PILGetReal("bandwidth", bandwidth))) {
    sprintf(msg, "Error reading the 'bandwidth' parameter!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Convert angle from [arc min] to [rad]
  *bandwidth = *bandwidth*M_PI/(60.*180.);

  return(status);
}




//////////////////////////
int generate_photons_main() 
{
  // Names of several input and output files:
  char orbit_filename[FILENAME_LENGTH];      // input: orbit
  char attitude_filename[FILENAME_LENGTH];   // input: attitude
  char photonlist_filename[FILENAME_LENGTH]; // output: photon list

  double t0;        // start time of the photon generation
  double timespan;  //  time span of the photon generation
  double bandwidth; // (half) width of the preselection band 
                    // along the path of the telescope axis [rad]


  // Catalog with orbit and attitude data over a particular timespan
  struct Telescope *sat_catalog=NULL;     
  // Number of entries in the orbit list ( <= orbit_nrows)
  long sat_nentries;                      

  gsl_rng *gsl_random_g=NULL; // pointer to GSL random number generator

  char msg[MAXMSG];           // error message buffer
  int status = EXIT_SUCCESS;  // error status flag


  // register HEATOOL
  set_toolname("generate_photons");
  set_toolversion("0.01");


  do {  // Beginning of ERROR HANDLING Loop.

    if ((status = generate_photons_getpar(orbit_filename, 
					  attitude_filename,
					  photonlist_filename, 
					  &t0, &timespan, &bandwidth))) 
      break;


    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);
    gsl_rng_env_setup();
    gsl_random_g = gsl_rng_alloc(gsl_rng_default);

    
    // Get the satellite catalog with the orbit and (telescope) attitude data:
    if ((status=get_satellite_catalog(&sat_catalog, &sat_nentries, t0, timespan, 
				      orbit_filename, attitude_filename))
	!=EXIT_SUCCESS) break;



  } while(0); // END of ERROR HANDLING Loop.



  // --- clean up ---

  // Release HEADAS random number generator:
  HDmtFree();
  gsl_rng_free(gsl_random_g);


  // Release memory of orbit/attitude catalog
  if (sat_catalog) free(sat_catalog);


  return(status);
}





