#include "measurement.h"


//#define EXTERN_PHOTON_LIST 1
void read_photon_list(struct Photon_Entry** pl, struct Detector);



////////////////////////////////////
// Main procedure.
int measurement_main() {
  double t0;        // starting time of the simulation
  double timespan;  // time span of the simulation
  double bandwidth; // (half) width of the preselection band 
                    // along the orbit [rad]

  char orbit_filename[FILENAME_LENGTH];   // filename of orbit file
  char attitude_filename[FILENAME_LENGTH];// filename of the attitude file
  long sat_nentries;                      // number of entries in the orbit array 
                                          // ( <= orbit_nrows)
  struct Telescope *sat_catalog=NULL;     // catalog with orbit and attitude data 
                                          // over a certain timespan
  struct vector n;                        // normalized vector perpendicular to 
                                          // the orbit plane

  int n_sourcefiles;                      // number of source input-files
  char source_filename[MAX_NSOURCEFILES][FILENAME_LENGTH]; // filenames of the 
                                          // individual source files
  fitsfile *source_catalog_files[MAX_NSOURCEFILES]; // file pointers to the 
                                          // individual source catalog files (FITS)
  int source_data_columns[5][3];          // column numbers of r.a., declination and
                                          // count rate in the individual files
  long nsources_pre;                      // number of totally available sources 
                                          // (ROSAT + RND + ...) 
                                          // in entire catalog and
                                          // in preselected catalog respectively

  long source_counter;                    // counter for source catalog access
  struct source_cat_entry *selected_catalog=NULL; // catalog of preselected sources
  struct source_cat_entry background;     // background properties: spectrum

  struct Photon_Entry *photon_list=NULL;  // photon list containing all actually 
                                          // created photons in the sky
  struct Photon_Entry *pl_entry=NULL;     // "counter" variable for the photon list

  struct Spectrum_Store spectrum_store;   // storage for different source spectra 
                                          // (including background spectrum)

  char spectrum_filename[N_SPECTRA_FILES][FILENAME_LENGTH]; 
  char rmf_name[FILENAME_LENGTH]; // FITS file containing the 
                                  //detector redistribution matrix (RMF)

  struct Telescope telescope;   // Telescope data (like FOV diameter or focal length)
  struct Detector detector;     // Detector data structure (containing the 
                                // pixel array, its width, ...)
  struct PSF_Store psf_store;   // Storage for the PSF (Point Spread Function) data 
                                // (for different off-axis angles and energies)
  char psf_filename[FILENAME_LENGTH]; // PSF input file

  struct Event_List_File event_list_file;

  gsl_rng *gsl_random_g;        // pointer to GSL random number generator
  int counter;                  // counter for several purpose

  char msg[MAXMSG];             // error output buffer
  int status=EXIT_SUCCESS;      // error status



  // register HEATOOL
  set_toolname("measurement");
  set_toolversion("0.01");

  // read parameters using PIL library
  status=measurement_getpar(orbit_filename, attitude_filename, &n_sourcefiles, 
			    source_filename, psf_filename, 
			    rmf_name, spectrum_filename[0], event_list_file.filename, 
			    &t0, &timespan, &telescope,
			    &detector, &bandwidth, &background.rate);


  // If an error has occurred during parameter input, break the program execution.
  if (status!=EXIT_SUCCESS) return(status);
  headas_chat(5,"\n");




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

  // typical timespan the satellite needs to slew over the entire FOV [s]:
  //  const double t_trans = 90.*60.*telescope.fov_diameter/(2.*M_PI);  

  // time step for sky scanning loop
  const double dt = 0.001;    // 10.* detector.integration_time; 




  // Initialize HEADAS random number generator and GSL generator for 
  // Gaussian distribution.
  HDmtInit(1);
  gsl_rng_env_setup();
  gsl_random_g = gsl_rng_alloc(gsl_rng_default);


  do {   // beginning of the error handling loop (will at most be run once)

    // get the satellite catalog with the orbit and (telescope) attitude data:
    if ((status=get_satellite_catalog(&sat_catalog, &sat_nentries, t0, timespan, 
				      orbit_filename, attitude_filename))
	!=EXIT_SUCCESS) break;



    // Set initial DETECTOR CONFIGURATION.

    // GENERAL SETTINGS
    // Calculate initial parameter values from the PIL parameters:
    detector.offset = detector.width/2;
    
    // Rescaling:
    detector.pixelwidth = detector.pixelwidth*1.e-6; // [mu m] -> [m]
    detector.ccsigma = detector.ccsigma*1.e-6;       // [mu m] -> [m]
    // size of the charge cloud [real pixel]
    detector.ccsize = 3.*detector.ccsigma/detector.pixelwidth;       
    // angle from [arc min] to [rad]
    telescope.fov_diameter = telescope.fov_diameter*M_PI/(60.*180.); 
    // convert angle from [arc min] to [rad]
    bandwidth = bandwidth*M_PI/(60.*180.);

    // Set the current detector frame to its initial value:
    detector.frame = -1;
    
    // DETECTOR SPECIFIC SETTINGS
    if (detector.type == FRAMESTORE) {
      headas_chat(5, "--> FRAMESTORE <--\n");

      // Set the first readout time such that the first readout is performed 
      // immediately at the beginning of the simulation (FRAMESTORE).
      detector.readout_time = t0 + detector.integration_time;
      
      detector.detector_action = framestore_detector_action;

      get_detector(&detector);
      
    } else if (detector.type == DEPFET) {
      headas_chat(5, "--> DEPFET <--\n");
      
      // Set the first readout time such that the first readout is performed 
      // immediately at the beginning of the simulation (DEPFET).
      detector.readout_time = t0 - detector.dead_time; 
      // The readout process starts at the center of the WFI detector, 
      // but for that purpose the current line has to be set to 0, so that the
      // first readout is performed in the middle of the detector array.
      detector.readout_line = 0;
      
      detector.detector_action = depfet_detector_action;
      
      headas_chat(5, "dead time: %lf\n", detector.dead_time);
      headas_chat(5, "clear time: %lf\n", detector.clear_time);
      headas_chat(5, "readout time: %lf\n", detector.readout_time);

      get_detector(&detector);
      
    } else if (detector.type == TES) {
      headas_chat(5, "--> TES Microcalorimeter Array <--\n");
    
      detector.detector_action = tes_detector_action;

      get_detector(&detector);

    } else if (detector.type == HTRS) {
      headas_chat(5, "--> HTRS <--\n");

      detector.detector_action = htrs_detector_action;

      htrs_get_detector(&detector);

    } else {

      status=EXIT_FAILURE;
      sprintf(msg, "Error: invalid detector type!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }
    

    // Consistency check for size of charge cloud:
    if (detector.ccsize > 1.) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: charge cloud size greater than pixel width!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }

    // END of DETECTOR CONFIGURATION SETUP



    // get the energy bins of the PHA channels
    if ((status=get_ebounds(&detector.ebounds, &detector.Nchannels, rmf_name))
	!=EXIT_SUCCESS) break;

    // get memory for the detector array
    //    if ((status=get_detector(&detector))!=EXIT_SUCCESS) break;


    // get the PSF
    if ((status=get_psf(&psf_store, psf_filename, &status))!=EXIT_SUCCESS) break;

    // TODO can be removed
    //    plot_array(psf_store.psf[0].data, psf_store.width, 1.0, "psf_image.png");

    // get the source spectra ( has to be called after get_psf() )
    if ((status=get_spectra(&spectrum_store, detector.Nchannels, 
			    spectrum_filename, N_SPECTRA_FILES, psf_store)) 
	!= EXIT_SUCCESS) break;

    // get the detector redistribution matrix (RMF)
    if ((status=get_rmf(&detector, rmf_name)) != EXIT_SUCCESS) break;


    // get the source catalogs
    if ((status=get_source_catalogs(&selected_catalog, n_sourcefiles, 
				    source_catalog_files, source_data_columns, 
				    source_filename))!=EXIT_SUCCESS) break;

    // create background data
    if ((background.rate*detector.integration_time*detector.width*detector.width<1.)
	&&(background.rate > 0.)) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: detector background count rate has to be larger "
	      "than %lf!\n", 
	      1./(detector.integration_time*detector.width*detector.width));
      HD_ERROR_THROW(msg,status);
      return(status);      
    }
    background.spectrum = &(spectrum_store.spectrum[0]);


    // Print some debug information:
    headas_chat(5, "source spectrum from '%s'\n", spectrum_filename[0]);
    headas_chat(5, "detector pixel width: %lf m\n", detector.pixelwidth);
    headas_chat(5, "charge cloud size: %lf pixel\n", detector.ccsize);
    headas_chat(5, "number of PHA channels: %d\n", detector.Nchannels);
    headas_chat(5, "lower PHA threshold: %ld\n\n", detector.low_threshold);



    // Create event list FITS file:
    headas_chat(5, "create FITS file '%s' for output of event list ...\n", 
		event_list_file.filename);
    // delete old event list
    remove(event_list_file.filename);

    // Create a new FITS file and a table for the event list.
    if (fits_create_file(&event_list_file.fptr, event_list_file.filename, &status)) 
      break;
    if (create_event_list_file(&event_list_file, detector, t0, t0+timespan, 

		        // HEADER keywords for event list FITS file:
			// telescope   CCD       instrument
			  "eROSITA",  "pnCCD1", "eROSITA",  &status)) break;



    //////////////////////

    // LOOP over all timesteps given the specified timespan from t0 to t0+timespan
    double time;                   // actual time
    //    double t_last_readout=t0;// time of the last detector readout
    double t_last=t0;              // last time of source preselection
    long sat_counter=0;            // counter for orbit readout loop
    long last_sat_counter=0;       // stores sat_counter of former repetition, 
                                   // so the searching loop
                                   // doesn't have to start at 0 every time.
    clear_detector(detector);      // clear the detector array (at the beginning 
                                   // there are no charges)

    // Beginning of actual simulation (after loading required data):
    headas_chat(5, "start measurement simulation ...\n");


    // Time loop:
    // Timesteps are typically a fraction (e.g. 1/10) of the time, the satellite 
    // takes to slew over the entire FOV.
    for(time=t0; (time<t0+timespan)&&(status==EXIT_SUCCESS); ) {

      // Print the current time step to STDOUT.
      //      headas_chat(5, "time: %lf\n", time);

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
      t_last= time;  // set the time of the last source preselection



      // PRESELECTION
      // Preselection of sources from the entire catalog to 
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


/*
      // LIGHTCURVES
      // Create lightcurves for new sources close to / within the FOV:

      // Determine telescope data (attitude):
      // Therefore we have to interpolate between the nearest data points known 
      // from the attitude file.
      telescope.nz = 
	normalize_vector(interpolate_vec(sat_catalog[sat_counter].nz, 
					 sat_catalog[sat_counter].time, 
					 sat_catalog[sat_counter+1].nz, 
					 sat_catalog[sat_counter+1].time, time));

      // Scan the preselected source catalog array for sources close to/within the FOV.
      for (source_counter=0; source_counter<nsources_pre; source_counter++) {
	// Compare each source to the unit vector specifiing the direction of 
	// the telescope:
	if (check_fov(selected_catalog[source_counter].r, telescope.nz, 
		      close_fov_min_align) == 0) {
	  // Source is close to the FOV, i.e. may be visible during the 
*/	  // integration timespan.
       /* // So check, whether a lightcurve was created recently for this source. 
          // If not, do this.
	  if (selected_catalog[source_counter]->lightcurve == NULL) {
	    // create lightcurve
	    if ((status=create_lightcurve(selected_catalog[source_counter], time, gsl_random_g)) != EXIT_SUCCESS) break;
	  } else if (selected_catalog[source_counter]->lightcurve[N_LIGHTCURVE_BINS].t < time) {
	    // create lightcurve
	    if ((status=create_lightcurve(selected_catalog[source_counter], time, gsl_random_g)) != EXIT_SUCCESS) break;
	  } */
/*	} else { // This part of code can be activated to save more memory 
	         // (The lightcurve free-memory is actually already implemented 
	         // in the preselection, 
	         // but here it can further reduce the amount of current used memory.)
	  // The source is not close to the FOV, so delete its lightcurve to save 
	  // memory (if one is available).
	  if (selected_catalog[source_counter].lightcurve != NULL) {
	    free(selected_catalog[source_counter].lightcurve);
	    selected_catalog[source_counter].lightcurve = NULL;
	    selected_catalog[source_counter].t_last_photon = -1.;
	  }
	}
      }  // END of loop over all sources in the close neighborhood of the 
         // telescope pointing direction.
      if (status != EXIT_SUCCESS) { break; }
      // END of lightcurve creation.
*/ // light curve creation is actually done in photon creation


      // CREATE PHOTONS for all sources close to the FOV
      for ( ; ((time<t_last+dt)&&(time<t0+timespan))&&(status==EXIT_SUCCESS); 
	    time+=dt /*detector.integration_time*/) {
	for (source_counter=0; source_counter<nsources_pre; source_counter++) {
	  // Check whether the source is inside the FOV:
	  // First determine telescope pointing direction at the actual time.
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
#ifndef EXTERN_PHOTON_LIST
	    if ((status=create_photons(&selected_catalog[source_counter], 
				       time, dt /*detector.integration_time*/, 
				       &photon_list, detector, gsl_random_g))
				       !=EXIT_SUCCESS) break; 
#else
	    if (photon_list == NULL) read_photon_list(&photon_list, detector);
#endif
	  }
	}
      } // End of fast inner time loop (photon creation and scanning of photon list).



      // SCAN PHOTON LIST
      // Perform the actual measurement (i.e., scan the photon list).
      sat_counter = last_sat_counter; 
      // Because we are now searching for events in the time interval [time-dt,time].
      while ((photon_list != NULL) && (status == EXIT_SUCCESS)) {
	// Call the detector action routine: this routine checks, whether the 
	// integration time is exceeded and performs the readout in this case. 
	// Otherwise it will simply do nothing.
        detector.detector_action(&detector, photon_list->photon.time, background, 
				 &event_list_file, &status);
	if (status != EXIT_SUCCESS) break;

	// If all photons up to the actual time have been treated, break the loop.
	if (photon_list->photon.time > time) {
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
	// First determine telescope pointing direction at the actual time.
	telescope.nz = 
	  normalize_vector(interpolate_vec(sat_catalog[sat_counter].nz, 
					   sat_catalog[sat_counter].time, 
					   sat_catalog[sat_counter+1].nz, 
					   sat_catalog[sat_counter+1].time, 
					   photon_list->photon.time));

	// Compare the photon direction to the unit vector specifiing the 
	// direction of the telescope axis:
#ifndef EXTERN_PHOTON_LIST
	if(check_fov(photon_list->photon.direction,telescope.nz,fov_min_align)==0) {
	  // Photon is inside the FOV.
#endif
	  
	  // Determine telescope data like direction etc. (attitude):
	  // Calculate nx: perpendicular to telescope axis and in the direction of
	  // the satellite motion:
	  telescope.nx = 
	    normalize_vector(interpolate_vec(sat_catalog[sat_counter].nx, 
					     sat_catalog[sat_counter].time, 
					     sat_catalog[sat_counter+1].nx, 
					     sat_catalog[sat_counter+1].time, 
					     photon_list->photon.time));
	  // Remove the component along the vertical direction nz 
	  // (nx must be perpendicular to nz!):
	  double scp = scalar_product(telescope.nz,telescope.nx);
	  telescope.nx.x -= scp*telescope.nz.x;
	  telescope.nx.y -= scp*telescope.nz.y;
	  telescope.nx.z -= scp*telescope.nz.z;
	  telescope.nx = normalize_vector(telescope.nx);

	  // ny is perpendicular to telescope axis and nx:
	  telescope.ny = normalize_vector(vector_product(telescope.nz, telescope.nx));
	  
	  // Measure the photon in a detector pixel, i.e., create the 
	  // corresponding charge there.
	  status =
	    measure(photon_list->photon, detector, telescope, psf_store, &event_list_file);
#ifndef EXTERN_PHOTON_LIST
	}
#endif

	// Move to the next entry in the photon list and clear the current entry.
	pl_entry = photon_list->next_entry;
	free(photon_list);
	photon_list = pl_entry;

      }  // END of scanning the photon list.

      last_sat_counter = sat_counter;
    }   // END of outer time loop.

  } while(0);  // END of the error handling loop.




  // --- cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // release HEADAS random number generator
  HDmtFree();
  gsl_rng_free(gsl_random_g);

  // Close the event list FITS file.
  fits_write_key(event_list_file.fptr, TLONG, "NPHOTONS", &photon_counter,
		 "number of measured photons", &status);
  if (event_list_file.fptr) fits_close_file(event_list_file.fptr, &status);

  // clear photon list
  clear_photon_list(&photon_list);

  // release memory of orbit catalog
  if (sat_catalog) free(sat_catalog);

  // release the light curves for the individual sources
  for (source_counter=0; source_counter<nsources_pre; source_counter++) {
    if (selected_catalog[source_counter].lightcurve != NULL) {
      free(selected_catalog[source_counter].lightcurve);
      selected_catalog[source_counter].lightcurve = NULL;
      selected_catalog[source_counter].t_last_photon = -1.;
    }
  }
  
  // release source catalogs
  free_source_catalogs(source_catalog_files,n_sourcefiles,&selected_catalog,&status);
  
  // release source spectra
  free_spectra(&spectrum_store, N_SPECTRA_FILES);

  // release memory of detector array
  if (detector.pixel) {
    for (counter = 0; counter < detector.width; counter++) {
      if (detector.pixel[counter]) {
	free(detector.pixel[counter]);
      }
    }
    free(detector.pixel);
  }

  // release memory of detector Redistribution Matrix
  free_rmf(detector.rmf);

  // release memory of detector EBOUNDS
  free_ebounds(detector.ebounds);

  // release memory of PSF:
  free_psf_store(psf_store);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished\n\n");

  return(status);
}




////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int measurement_getpar(
		       char orbit_filename[],   // filename of the orbit file (FITS)
		       // filename of the attitude file (FITS)
		       char attitude_filename[],
		       int *n_sourcefiles,      // number of input source files
		       // array containing the filename of each source file
		       char source_filename[MAX_NSOURCEFILES][FILENAME_LENGTH], 
		       char psf_filename[],     // PSF FITS file
		       // FITS file containing the detector EBOUNDS
		       char rmf_name[],       
		       // PHA file containing the default source spectrum:
		       char spectrum_filename[],
		       // filename of the output evenlist (FITS file):
		       char eventlist_filename[],
		       double *t0,              // start time for the simulation
		       double *timespan,        // time span for the simulation
		       struct Telescope *telescope,
		       struct Detector *detector,
		   // width of the band around the sky for source preselection [rad]:
		       double *bandwidth,       
		       float *background_countrate
		       )
{
  char cbuffer[FILENAME_LENGTH];// buffer to access the different source filenames

  char msg[MAXMSG];             // error output buffer
  int status=EXIT_SUCCESS;      // error status


  // get the filename of the orbit file (FITS file)
  if ((status = PILGetFname("orbitfile", orbit_filename))) {
    sprintf(msg, "Error reading the filename of the orbitfile!\n");
    HD_ERROR_THROW(msg,status);
  }

  // get the filename of the attitude file (FITS file)
  else if ((status = PILGetFname("attitudefile", attitude_filename))) {
    sprintf(msg, "Error reading the filename of the attitudefile!\n");
    HD_ERROR_THROW(msg,status);
  }

  // get the number of source input-files
  else if ((status = PILGetInt("n_sourcefiles", n_sourcefiles))) {
    sprintf(msg, "Error reading the 'number of sourcefiles' parameter!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the filenames of the individual source catalogs.
  else {
    int input_counter;
    for(input_counter=0; input_counter<*n_sourcefiles; input_counter++) {
      sprintf(cbuffer,"sourcefile%d",input_counter+1);
      if ((status = PILGetFname(cbuffer, source_filename[input_counter]))) {
	sprintf(msg, "Error reading the name of the sourcefile No %d!\n", 
		input_counter);
	HD_ERROR_THROW(msg,status);
      }
    }
    if (status) return(status);
  }


  // Get the detector type (FRAMESTORE, DEPFET, TES microcalorimeter, ...)
  int type;
  if ((status = PILGetInt("detectortype", &type))) {
    sprintf(msg, "Error reading the detector type!\n");
    HD_ERROR_THROW(msg,status);
    return (status);
  } else {
    switch (type) {
      // According to the different detector types, the program
      // has to read different parameters.
    case 1: 
      detector->type = FRAMESTORE; 

      // Get the integration time for the FRAMESTORE CCD.
      if ((status = PILGetReal("integration_time", &detector->integration_time))) {
	sprintf(msg, "Error reading the integration time!\n");
	HD_ERROR_THROW(msg,status);
      }

      break;

    case 2: 
      detector->type = DEPFET;     

      // Get the dead time for the DEPFET APS (readout time per line).
      if ((status = PILGetReal("dead_time", &detector->dead_time))) {
	sprintf(msg, "Error reading the dead time!\n");
	HD_ERROR_THROW(msg,status);
      } 
      // Get the clear time for the DEPFET APS (time required to clear one line).
      else if ((status = PILGetReal("clear_time", &detector->clear_time))) {
	sprintf(msg, "Error reading the clear time!\n");
	HD_ERROR_THROW(msg,status);
      }
      
      break;

    case 3: 
      detector->type = TES;        
      
      break;

    default:     
      sprintf(msg, "Error: incorrect detector type!\n");
      HD_ERROR_THROW(msg,status);
    }
  }
  if (status) return(status);
  // END of handling different detector types.



  // get the filename of the PSF data file (FITS file)
  if ((status = PILGetFname("psffile", psf_filename))) {
    sprintf(msg, "Error reading the filename of the PSF file!\n");
    HD_ERROR_THROW(msg,status);
  }

  // filename of the detector redistribution file (FITS file)
  else if ((status = PILGetFname("rmf", rmf_name))) {
    sprintf(msg, "Error reading the filename of the detector" 
	    "redistribution matrix file (RMF)!\n");
    HD_ERROR_THROW(msg,status);
  }

  // filename of the default source spectrum (PHA FITS file)
  else if ((status = PILGetFname("spectrum", spectrum_filename))) {
    sprintf(msg, "Error reading the filename of the default spectrum (PHA)!\n");
    HD_ERROR_THROW(msg,status);
  }

  // filename of the output event list (FITS file)
  else if ((status = PILGetFname("eventlist", eventlist_filename))) {
    sprintf(msg, "Error reading the filename of the eventlist!\n");
    HD_ERROR_THROW(msg,status);
  }

  // get the start time of the simulation
  else if ((status = PILGetReal("t0", t0))) {
    sprintf(msg, "Error reading the 't0' parameter!\n");
    HD_ERROR_THROW(msg,status);
  }

  // get the timespan for the simulation
  else if ((status = PILGetReal("timespan", timespan))) {
    sprintf(msg, "Error reading the 'timespan' parameter!\n");
    HD_ERROR_THROW(msg,status);
  }

  // read the diameter of the FOV (in arcmin)
  else if ((status = PILGetReal("fov_diameter", &telescope->fov_diameter))) {
    sprintf(msg, "Error reading the diameter of the FOV!\n");
    HD_ERROR_THROW(msg,status);
  }

  // read the focal length [m]
  else if ((status = PILGetReal("focal_length", &telescope->focal_length))) {
    sprintf(msg, "Error reading the focal length!\n");
    HD_ERROR_THROW(msg,status);
  }

  // [pixel]
  else if ((status = PILGetInt("det_width", &detector->width))) {
    sprintf(msg, "Error reading the width of the detector!\n");
    HD_ERROR_THROW(msg,status);
  }

  // [mu m]
  else if ((status = PILGetReal("det_pixelwidth", &detector->pixelwidth))) {
    sprintf(msg, "Error reading the width of the detector pixels!\n");
    HD_ERROR_THROW(msg,status);
  }

  // [mu m]
  else if ((status = PILGetReal("ccsigma", &detector->ccsigma))) {
    sprintf(msg, "Error reading the charge cloud sigma!\n");
    HD_ERROR_THROW(msg,status);
  }
    
  // read charge cloud sigma
  else if ((status = PILGetReal4("background_countrate", background_countrate))) {
    sprintf(msg, "Error reading the background countrate!\n");
    HD_ERROR_THROW(msg,status);
  }

  // read the (half) width of the preselection band (in arcmin)
  else if ((status = PILGetReal("bandwidth", bandwidth))) {
    sprintf(msg, "Error reading the 'bandwidth' parameter!\n");
    HD_ERROR_THROW(msg,status);
  }

  return(status);
}





//////////////////////////////////////////////////////////////////////////////
// Function performs a measurement for a specific source inside dthe FOV 
// with given telescope data (direction of axis, direction of motion).
int measure(
	    struct Photon photon, 
	    struct Detector detector,
	    struct Telescope telescope,
	    struct PSF_Store psf_store,
	    struct Event_List_File* event_list_file
	    ) 
{
#ifndef EXTERN_PHOTON_LIST
  struct Point2d position;  // Photon impact position on the detector in [m]
  
  int status=EXIT_SUCCESS;


  // Convolution with PSF: get the position [m], where the photon hits the detector.
  // Function returns 0, if the photon does not fall on the detector. 
  // If it hits a detector pixel, the return value is 1.
  if (get_psf_pos(&position, photon, telescope, psf_store)) {
    // Check whether the photon hits the detector within the FOV. (By the effects 
    // of the mirrors it might have been scattered over the edge of the FOV, 
    // although the source is inside the FOV.)
    if (sqrt(pow(position.x,2.)+pow(position.y,2.)) < 
	tan(telescope.fov_diameter)*telescope.focal_length) {

      // Increase the counter of measured photons:
      photon_counter++;

      int x[4], y[4];
      double fraction[4];
      int npixels = get_pixel_square(detector, position, x, y, fraction);
      
      if ((detector.type == FRAMESTORE) || (detector.type == DEPFET)) {
	// Determine a detector channel (PHA channel) according to RMF.
	long channel = detector_rmf(photon.energy, detector.rmf);
	// Get the corresponding created charge.
	float charge = get_charge(channel, detector.ebounds);

	if (channel <= 0) printf("ERROR CHANNEL <= 0!\n");  // TODO
	if (channel > 4096) printf("ERROR CHANNEL too large!\n");

	int count;
	for (count=0; count<npixels; count++) {
	  if (x[count] != INVALID_PIXEL) {
	    detector.pixel[x[count]][y[count]].charge += 
	      charge * fraction[count] * 
	      // |        |-> charge fraction due to split events
	      // |-> charge created by incident photon
	      detector_active(x[count], y[count], detector, photon.time);
	    // |-> "1" if pixel can measure charge, "0" else
	  }
	}

      } else if (detector.type == TES) {
	struct Event event;

	// Store the photon charge and the new arrival time:
	event.pha = get_pha(photon.energy, detector);  // TODO: RMF
	event.time = photon.time;
	event.xi = x[0];
	event.yi = y[0];
	event.grade = 0;
	event.frame = detector.frame;

	// Add the event to the FITS event list.
	if (event.pha >= detector.low_threshold) { // Check lower PHA threshold
	  // There is an event in this pixel, so insert it into eventlist:
	  add_eventtbl_row(event_list_file, event, &status);
	}

      } // END of detector type TES

    } // END of Check whether hitting point is inside the FOV
  } // END of 'get_psf_pos(...)'

  return(status);

#else
  double det_x, det_y;   // Photon impact position on the detector in [m]
  
  det_x = photon.direction.x - 8*75.e-6;
  det_y = photon.direction.y - 8*75.e-6;
  float charge = photon.energy;

  // Increase the counter of measured photons:
  photon_counter++;

  // Calculate integer detector pixels: 
  // It is important to perform the (int) operation after adding the offset!
  // Otherwise there may be errors for negativ 'det_x' or 'det_y' values.
  int det_xi = (int)(det_x/detector.pixelwidth + (double)(detector.width/2));
  int det_yi = (int)(det_y/detector.pixelwidth + (double)(detector.width/2));

  // Calculate split events due to the spread of the charge cloud:
  double partition[3][3];
  split_events(det_x, det_y, det_xi, det_yi, detector, partition); 

  // Add the charge which is created by the photon to the 
  // detector pixels considering split events.
  int countx, county;
  const int countx_max=MAX(3,4+det_xi-detector.width);
  const int county_max=MAX(3,4+det_yi-detector.width);

  for (countx=MAX(0,1-det_xi); countx<countx_max; countx++) {
    for (county=MAX(0,1-det_yi); county<county_max; county++) {
      if(((det_xi+countx-1)<detector.width)&&((det_yi+county-1)<detector.width)) {
	detector.pixel[det_xi+countx-1][det_yi+county-1].charge += 
	  charge * partition[countx][county] * 
	  // |        |-> charge fraction due to split events
	  // |-> charge created by incident photon
	  detector_active(det_xi+countx-1, det_yi+county-1, detector, photon.time);
	// |-> "1" if pixel can measure charge, "0" else
      }
    }
  }
  
#endif
}




void read_photon_list(struct Photon_Entry** pl, struct Detector detector)
{
  fitsfile *fptr;
  long nrows, row;
  char filename[] = "Createdphotonlist_window_comp_relInt_00000000.010000.fits";

  int status = EXIT_SUCCESS;
  
  do {
    headas_chat(5, "open FITS file '%s' containing the photon list ...\n", filename);
    if (fits_open_table(&fptr, filename, READONLY, &status)) break;

    if (fits_get_num_rows(fptr, &nrows, &status)) break;

    double pha_buffer;
    long pha;
    double x, y;
    double time;
    int anynul = 0;
    for (row = 0; row < nrows; row++) {
      fits_read_col(fptr, TDOUBLE, 1, row+1,1,1, &pha_buffer, &pha_buffer, 
		    &anynul, &status);
      pha = (long)pha_buffer;
      fits_read_col(fptr, TDOUBLE, 2, row+1,1,1, &time, &time, &anynul, &status);
      fits_read_col(fptr, TDOUBLE, 5, row+1,1,1, &x, &x, &anynul, &status);
      fits_read_col(fptr, TDOUBLE, 6, row+1,1,1, &y, &y, &anynul, &status);

      struct Photon new_photon;
      new_photon.direction.x = x;
      new_photon.direction.y = y;
      new_photon.energy = get_charge(pha, detector.ebounds);
      new_photon.time = time;

      insert_photon(pl, new_photon);
    }

  } while (0);

  fits_close_file(fptr, &status);

  if (status != EXIT_SUCCESS) printf("Error while reading the photon list!\n");
}


