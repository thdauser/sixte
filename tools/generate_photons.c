#include "generate_photons.h"



////////////////////////////
int generate_photons_getpar(
			    char orbit_filename[],
			    char attitude_filename[],
			    // number of input source catalog files
			    int *n_sourcefiles,  
			    // array containing the filename of each source file
			    char source_filename[MAX_NSOURCEFILES][FILENAME_LENGTH],
			    // PHA file containing the default source spectrum
			    char spectrum_filename[],
			    char rmf_filename[],
			    char photonlist_filename[],
			    double *t0,
			    double *timespan,
			    double *bandwidth,
			    struct Telescope *telescope
			    )
{
  // filename-buffer to access the different source files
  char cbuffer[FILENAME_LENGTH];

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

  // Get the number of source input-files
  else if ((status = PILGetInt("n_sourcefiles", n_sourcefiles))) {
    sprintf(msg, "Error reading the number of source catalog files!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the filenames of the individual source catalogs.
  else {
    int input_counter;
    for(input_counter=0; input_counter<MAX_NSOURCEFILES; input_counter++) {

      sprintf(cbuffer,"sourcefile%d",input_counter+1);

      if (input_counter<*n_sourcefiles) {
	// Read the source file from using PIL.
	if ((status = PILGetFname(cbuffer, source_filename[input_counter]))) {
	  sprintf(msg, "Error reading the name of the sourcefile No %d!\n", 
		  input_counter);
	  HD_ERROR_THROW(msg,status);
	}

      } else {
	// Fill redundant input slots for source files with NULL value,
	// in order to avoid errors with the HD_PARSTAMP routine.
	PILPutFname(cbuffer, "");

      }
	
    }
    if (status) return(status);
  }

  // Get the filename of the default source spectrum (PHA FITS file)
  if ((status = PILGetFname("spectrum", spectrum_filename))) {
    sprintf(msg, "Error reading the filename of the default spectrum (PHA)!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the filename of the detector redistribution file (FITS file)
  else if ((status = PILGetFname("rmf", rmf_filename))) {
    sprintf(msg, "Error reading the filename of the detector" 
	    "redistribution matrix file (RMF)!\n");
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

  // Get the diameter of the FOV (in arcmin)
  else if ((status = PILGetReal("fov_diameter", &telescope->fov_diameter))) {
    sprintf(msg, "Error reading the diameter of the FOV!\n");
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
  char spectrum_filename[N_SPECTRA_FILES][FILENAME_LENGTH];// input: source spectra
  char rmf_filename[FILENAME_LENGTH];        // input: detector RMF
  char photonlist_filename[FILENAME_LENGTH]; // output: photon list
  
  // Several input source catalog files:
  int n_sourcefiles;   // number of input source files
  // Filenames of the individual source catalog files (FITS)
  char source_filename[MAX_NSOURCEFILES][FILENAME_LENGTH]; 

  double t0;        // start time of the photon generation
  double timespan;  //  time span of the photon generation
  double bandwidth; // (half) width of the preselection band 
                    // along the path of the telescope axis [rad]

  // Catalog with orbit and attitude data over a particular timespan
  struct Telescope *sat_catalog=NULL;     
  // Number of entries in the orbit list ( <= orbit_nrows)
  long sat_nentries;                      

  // Source catalogs (FITS files)
  fitsfile *source_catalog_files[MAX_NSOURCEFILES];
  // Catalog of preselected sources along the path of the telescope axis
  struct source_cat_entry *selected_catalog=NULL;
  // Column numbers of r.a., declination and count rate in the individual files
  int source_data_columns[5][3];       
  // Number of totally available sources (ROSAT + RND + ...) in entire catalog and
  // in preselected catalog respectively.
  long source_counter, nsources_pre=0;

  // Storage for different source spectra (including background spectrum).
  struct Spectrum_Store spectrum_store; 

  // Telescope data (like FOV diameter or focal length)
  struct Telescope telescope; 
  // Detector data structure (containing the pixel array, its width, 
  // RMF, EBOUNDS ...)
  struct Detector detector;   

  gsl_rng *gsl_random_g=NULL; // pointer to GSL random number generator

  char msg[MAXMSG];           // error message buffer
  int status = EXIT_SUCCESS;  // error status flag


  // register HEATOOL
  set_toolname("generate_photons");
  set_toolversion("0.01");


  do {  // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    if ((status = generate_photons_getpar(orbit_filename, attitude_filename,
					  &n_sourcefiles, source_filename,
					  spectrum_filename[0], rmf_filename, 
					  photonlist_filename,
					  &t0, &timespan, &bandwidth,
					  &telescope))) break;


    // Set last_update to such a small value, that a preselection of the 
    // source catalog is performed at the first timestep (last_update 
    // contains the time of the last source catalog preselection.):
    const double former_time = t0 - ORBIT_UPDATE_TIME - 100.;
    double last_update = former_time;

    // Defines the mathematical meaning of 'close' in the context that for 
    // sources 'close to the FOV' the simulation creates a light curve.
    const double close_mult = 1.5; 
  
    // Calculate the minimum cos-value for sources inside the FOV and 
    // in the preselection band respectively:
    
    // minimum cos-value for sources inside the FOV angle(x0,source) <= 1/2 * diameter
    const double fov_min_align = cos(telescope.fov_diameter/2.); 
    
    // minimum cos-value for sources close to the FOV (in the direct neighborhood) 
    const double close_fov_min_align = cos(close_mult*telescope.fov_diameter/2.); 

    // maximum cos-value (minimum sin-value) for sources within the 
    // preselection band along the orbit (angle(n,source) > 90-bandwidth)
    const double pre_max_align = sin(bandwidth);

    // time step for sky scanning loop
    const double dt = 0.001;    // 10.* detector.integration_time; 


    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);
    gsl_rng_env_setup();
    gsl_random_g = gsl_rng_alloc(gsl_rng_default);

    
    // Get the satellite catalog with the orbit and (telescope) attitude data:
    if ((status=get_satellite_catalog(&sat_catalog, &sat_nentries, t0, 
				      timespan, orbit_filename, 
				      attitude_filename)) !=EXIT_SUCCESS) break;

    // Get the energy bins of the PHA channels:
    if ((status=get_ebounds(&detector.ebounds, &detector.Nchannels, rmf_filename))
	!=EXIT_SUCCESS) break;

    // Get the source spectra:
    if ((status=get_spectra(&spectrum_store, detector.Nchannels, spectrum_filename,
			    N_SPECTRA_FILES)) != EXIT_SUCCESS) break;
        
    // Get the source catalogs:
    if ((status=get_source_catalogs(&selected_catalog, n_sourcefiles, 
				    source_catalog_files, source_data_columns, 
				    source_filename))!=EXIT_SUCCESS) break;


    // ---- End of Initialization ----

  } while(0); // END of ERROR HANDLING Loop.



  // --- Clean up ---

  // Release HEADAS random number generator:
  HDmtFree();
  gsl_rng_free(gsl_random_g);


  // Release memory of orbit/attitude catalog
  if (sat_catalog) free(sat_catalog);

  // Release the light curves for the individual sources
  for (source_counter=0; source_counter<nsources_pre; source_counter++) {
    if (selected_catalog[source_counter].lightcurve != NULL) {
      free(selected_catalog[source_counter].lightcurve);
      selected_catalog[source_counter].lightcurve = NULL;
      selected_catalog[source_counter].t_last_photon = -1.;
    }
  }
  
  // Release source catalogs
  free_source_catalogs(source_catalog_files,n_sourcefiles,
		       &selected_catalog,&status);
  
  // Release source spectra
  free_spectra(&spectrum_store, N_SPECTRA_FILES);

  return(status);
}





