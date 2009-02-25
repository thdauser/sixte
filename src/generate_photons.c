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
    HD_ERROR_THROW(msg, status);
  }

  // Get the filename of the Attitude file (FITS file):
  else if ((status = PILGetFname("attitude_filename", attitude_filename))) {
    sprintf(msg, "Error reading the filename of the attitude file!\n");
    HD_ERROR_THROW(msg, status);
  }

  // Get the number of source input-files
  else if ((status = PILGetInt("n_sourcefiles", n_sourcefiles))) {
    sprintf(msg, "Error reading the number of source catalog files!\n");
    HD_ERROR_THROW(msg, status);
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

  // Photon list containing all created photons in the sky
  struct Photon_Entry *photon_list=NULL;  

  // Pointer to the FITS file for the output for the photon list.
  fitsfile *photonlist_fptr = NULL;

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
    
    // minimum cos-value for sources inside the FOV 
    // angle(x0,source) <= 1/2 * diameter
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


    // Delete old photon list FITS file:
    remove(photonlist_filename);
    // Create a new FITS file for the output of the photon list:
    if ((create_photonlist_file(&photonlist_fptr, photonlist_filename, &status))) break;


    // --- End of Initialization ---



    // --- Beginning of Photon generation process ---

    // LOOP over all timesteps given the specified timespan from t0 to t0+timespan
    double time;             // current time

    long sat_counter=0;      // counter for orbit readout loop
    long last_sat_counter=0; // stores sat_counter of former repetition, 
                             // so the searching loop
                             // doesn't have to start at 0 every time.
    long photon_row=0;       // current row in photon list FITS file;

    struct vector n;   // normalized vector perpendicular to the orbital plane

    struct Photon_Entry *pl_entry=NULL;     // "counter" variable for the photon list


    // Beginning of actual simulation (after loading required data):
    headas_chat(5, "start photon generation process ...\n");


    // Time loop:
    // Timesteps are typically a fraction (e.g. 1/10) of the time, the satellite 
    // takes to slew over the entire FOV.
    for(time=t0; (time<t0+timespan)&&(status==EXIT_SUCCESS); time+=dt) {

      // Print the current time step to STDOUT.
      // headas_chat(5, "time: %lf\n", time);

      // Get the last orbit entry before the time 'time':
      // (in order to interpolate the actual position and velocity between 
      // the neighboring calculated orbit positions)
      for(; sat_counter<sat_nentries; sat_counter++) {
	if(sat_catalog[sat_counter].time>time) {
	  break;
	}
      }
      if(fabs(sat_catalog[--sat_counter].time-time)>600.) { 
	// no entry within 10 minutes !!
	status = EXIT_FAILURE;
	sprintf(msg, "Error: no adequate orbit entry for time %lf!\n", time);
	HD_ERROR_THROW(msg,status);
	break;
      }



      // PRESELECTION
      // Preselection of sources from the comprehensive catalog to 
      // improve the performance of the simulation:
      if (sat_catalog[sat_counter].time-last_update > ORBIT_UPDATE_TIME) {
	// Preselect sources from the entire source catalog according to the 
	// satellite's direction of movement.
	// Calculate normalized vector perpendicular to the orbit plane:
	n = 
	  normalize_vector(vector_product(normalize_vector(sat_catalog[sat_counter].r),
					  normalize_vector(sat_catalog[sat_counter].v)));
	if ((status=get_preselected_catalog(selected_catalog, &nsources_pre, 
					    n_sourcefiles, source_catalog_files, 
					    source_data_columns, n, pre_max_align, 
					    spectrum_store, N_SPECTRA_FILES))
	    !=EXIT_SUCCESS) break;

	// update the catalog-update-counter
	last_update = sat_catalog[sat_counter].time;
      }
      // END of preselection



      // CREATE PHOTONS for all sources  CLOSE TO  the FOV
      for (source_counter=0; source_counter<nsources_pre; source_counter++) {
	// Check whether the source is inside the FOV:
	
	// First determine telescope pointing direction at the actual time.
	// TODO: replace this calculation by orbit interpolation.
	telescope.nz =
	  normalize_vector(interpolate_vec(sat_catalog[sat_counter].nz, 
					   sat_catalog[sat_counter].time, 
					   sat_catalog[sat_counter+1].nz,
					   sat_catalog[sat_counter+1].time,time));
	
	
	// Compare the source direction to the unit vector specifiing the 
	// direction of the telescope:
	if(check_fov(selected_catalog[source_counter].r, telescope.nz, 
		     close_fov_min_align) == 0) {
	  // The source is inside the FOV  => create photons:
	  if ((status=create_photons(&selected_catalog[source_counter], time, dt, 
				     &photon_list, detector, gsl_random_g))
	      !=EXIT_SUCCESS) break; 
	}
      }



      // SCAN PHOTON LIST
      // Perform the actual measurement (i.e., scan the photon list).
      sat_counter = last_sat_counter; 
      // Because we are now searching for events in the time interval [time-dt,time].
      while ((photon_list != NULL) && (status == EXIT_SUCCESS)) {
	// If all photons up to the actual time have been treated, break the loop.
	if ((photon_list->photon.time > time + dt)||
	    (photon_list->photon.time > t0+timespan)) {
	  break;
	}

	// Get the last orbit entry before the time 'photon_list->photon.time'
	// (in order to interpolate the position and velocity at this time  between 
	// the neighboring calculated orbit positions):
	for( ; sat_counter<sat_nentries; sat_counter++) {
	  if(sat_catalog[sat_counter].time>photon_list->photon.time) {
	    break;
	  }
	}
	if(fabs(sat_catalog[--sat_counter].time-time)>600.) { 
	  // no entry within 10 minutes !!
	  status = EXIT_FAILURE;
	  sprintf(msg, "Error: no adequate orbit entry for time %lf!\n", time);
	  HD_ERROR_THROW(msg,status);
	  break;
	}

	// Check whether the photon is inside the FOV:
	// First determine telescope pointing direction at the current time.
	telescope.nz = 
	  normalize_vector(interpolate_vec(sat_catalog[sat_counter].nz, 
					   sat_catalog[sat_counter].time, 
					   sat_catalog[sat_counter+1].nz, 
					   sat_catalog[sat_counter+1].time, 
					   photon_list->photon.time));

	// Compare the photon direction to the unit vector specifiing the 
	// direction of the telescope axis:
	if(check_fov(photon_list->photon.direction,telescope.nz,fov_min_align)==0) {
	  // Photon is inside the FOV!

	  // Add the photon to the photon list file:
	  fits_insert_rows(photonlist_fptr, photon_row++, 1, &status);
	  fits_write_col(photonlist_fptr, TDOUBLE, 1, photon_row, 1, 1, 
			 &photon_list->photon.time, &status);
	  fits_write_col(photonlist_fptr, TFLOAT, 2, photon_row, 1, 1, 
			 &photon_list->photon.energy, &status);
	  fits_write_col(photonlist_fptr, TDOUBLE, 3, photon_row, 1, 1, 
			 &photon_list->photon.ra, &status);
	  fits_write_col(photonlist_fptr, TDOUBLE, 4, photon_row, 1, 1, 
			 &photon_list->photon.dec, &status);
	}

	// Move to the next entry in the photon list and clear the current entry.
	pl_entry = photon_list->next_entry;
	free(photon_list);
	photon_list = pl_entry;

      }  // END of scanning the photon list.

      last_sat_counter = sat_counter;
    }   // END of outer time loop.

    // --- End of photon creation process ---

  } while(0); // END of ERROR HANDLING Loop.



  // --- Clean up ---

  // Close FITS file
  if (photonlist_fptr) fits_close_file(photonlist_fptr, &status);

  // Release HEADAS random number generator:
  HDmtFree();
  gsl_rng_free(gsl_random_g);

  // clear photon list
  clear_photon_list(&photon_list);

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

  // release memory of detector EBOUNDS
  free_ebounds(detector.ebounds);


  return(status);
}





