#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include "photon_generation.h"
#include "clusters.c"


////////////////////////////
int photon_generation_getpar(
			     char orbit_filename[],
			     char attitude_filename[],
			     // number of input source catalog files
			     int* n_sourcefiles,  
			     // array containing the filename of each source file
			     char** source_filename,
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
    for(input_counter=0; input_counter<MAX_N_POINTSOURCEFILES; input_counter++) {

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


  // Convert angles from [arc min] to [rad]
  *bandwidth = *bandwidth*M_PI/(60.*180.);
  telescope->fov_diameter = telescope->fov_diameter*M_PI/(180.*60.);

  return(status);
}




//////////////////////////
int photon_generation_main() 
{
  // Names of several input and output files:
  char orbit_filename[FILENAME_LENGTH];      // input: orbit
  char attitude_filename[FILENAME_LENGTH];   // input: attitude
  char spectrum_filename[N_SPECTRA_FILES][FILENAME_LENGTH]; // input: source spectra
  char rmf_filename[FILENAME_LENGTH];        // input: detector RMF
  char photonlist_filename[FILENAME_LENGTH]; // output: photon list
  
  // Several input source catalog files:
  int n_sourcefiles;   // number of input source files
  // Filenames of the individual source catalog files (FITS):
  char** source_filename=NULL;
  // X-ray Cluster image:
  char cluster_filename[FILENAME_LENGTH];    // input: cluster image file
  ClusterImage* cluster_image=NULL;

  // New data structures for point sources:
  PointSourceFiles* pointsourcefiles=NULL;
  PointSourceCatalog* pointsourcecatalog=NULL;
  long source_counter;

  double t0;        // start time of the photon generation
  double timespan;  // time  span of the photon generation
  double bandwidth; // (half) width of the preselection band 
                    // along the path of the telescope axis [rad]

  // Catalog with orbit and attitude data over a particular timespan
  struct Telescope *sat_catalog=NULL;     
  // Number of entries in the orbit list ( <= orbit_nrows)
  long sat_nentries;                      

  // Storage for different source spectra (including background spectrum).
  struct Spectrum_Store spectrum_store; 

  // Telescope data (like FOV diameter or focal length)
  struct Telescope telescope; 
  // Detector data structure (containing the pixel array, its width, 
  // RMF, EBOUNDS ...)
  Detector* detector=NULL;

  // Photon list containing all created photons in the sky
  struct Photon_Entry *photon_list=NULL;  

  // Pointer to the FITS file for the output for the photon list.
  fitsfile *photonlist_fptr = NULL;

  gsl_rng *gsl_random_g=NULL; // pointer to GSL random number generator

  char msg[MAXMSG];           // error message buffer
  int status = EXIT_SUCCESS;  // error status flag


  // register HEATOOL
  set_toolname("photon_generation");
  set_toolversion("0.01");


  do {  // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    detector = get_Detector(&status);
    if(status!=EXIT_SUCCESS) break;
    
    int count;
    source_filename=(char**)malloc(MAX_N_POINTSOURCEFILES*sizeof(char*));
    if (source_filename!=NULL) {
      for(count=0; (count<MAX_N_POINTSOURCEFILES)&&(status==EXIT_SUCCESS); count++) {
	source_filename[count] = (char*)malloc(FILENAME_LENGTH*sizeof(char));
	if(source_filename[count]==NULL) {
	  status=EXIT_FAILURE;
	  sprintf(msg, "Error: not enough memory!\n");
	  HD_ERROR_THROW(msg,status);
	  break;
	}
      }
    } else {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }


    if((status=photon_generation_getpar(orbit_filename, attitude_filename,
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
    const double dt = 0.1;


    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);
    gsl_rng_env_setup();
    gsl_random_g = gsl_rng_alloc(gsl_rng_default);

    
    // Get the satellite catalog with the orbit and (telescope) attitude data:
    if ((status=get_satellite_catalog(&sat_catalog, &sat_nentries, t0, 
				      timespan, orbit_filename, 
				      attitude_filename)) !=EXIT_SUCCESS) break;

    // Read the detector RMF and EBOUNDS from the specified file and assign them to the 
    // Detector data structure.
    if ((status=detector_assign_rsp(detector, rmf_filename)) != EXIT_SUCCESS) break;

    // Get the source spectra:
    if ((status=get_spectra(&spectrum_store, detector->rmf->NumberChannels, 
			    spectrum_filename, N_SPECTRA_FILES)) != EXIT_SUCCESS) break;
    
    // Get the source catalogs:
    pointsourcefiles = get_PointSourceFiles(n_sourcefiles, source_filename, &status);
    if (status != EXIT_SUCCESS) break;

    // Get the specified galaxy cluster image:
    if ((status = PILGetFname("cluster_filename", cluster_filename))) {
      sprintf(msg, "Error reading the filename of the cluster image file!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }
    if(strlen(cluster_filename) > 0) {
      cluster_image = get_ClusterImage_fromFile(cluster_filename, &status);
      if (status != EXIT_SUCCESS) break;
    }

    // Delete old photon list FITS file:
    remove(photonlist_filename);
    // Create a new FITS file for the output of the photon list:
    if ((create_photonlist_file(&photonlist_fptr, photonlist_filename, &status))) 
      break;


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

    struct Photon_Entry *pl_entry=NULL; // "counter" variable for the photon list


    // Beginning of actual simulation (after loading required data):
    headas_chat(5, "start photon generation process ...\n");


    // Time loop:
    // Timesteps are typically a fraction (e.g. 1/10) of the time, the satellite 
    // takes to slew over the entire FOV.
    for(time=t0; (time<t0+timespan)&&(status==EXIT_SUCCESS); time+=dt) {
      printf("time: %lf\n", time);

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
	if((status=get_PointSourceCatalog(pointsourcefiles, &pointsourcecatalog, n, 
					  pre_max_align, spectrum_store))
	   !=EXIT_SUCCESS) break;

	// Update the catalog-update-counter
	last_update = sat_catalog[sat_counter].time;
      }
      // END of preselection



      // CREATE PHOTONS for all POINT sources CLOSE TO the FOV
      for (source_counter=0; source_counter<pointsourcecatalog->nsources; 
	   source_counter++) {
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
	//if(check_fov(selected_catalog[source_counter].r, telescope.nz, 
	//	     close_fov_min_align) == 0) {
	struct vector source_vector = 
	  unit_vector(pointsourcecatalog->sources[source_counter].ra, 
		      pointsourcecatalog->sources[source_counter].dec);
	if (check_fov(&source_vector, &telescope.nz, close_fov_min_align) == 0) {

	  // The source is inside the FOV  => create photons:
	  //	  if ((status=create_photons(&selected_catalog[source_counter], time, dt, 
	  //				     &photon_list, detector, gsl_random_g))
	  //	      !=EXIT_SUCCESS) break; 
	  if ((status=create_photons(&(pointsourcecatalog->sources[source_counter]), 
				     time, dt, &photon_list, detector, gsl_random_g))
	      !=EXIT_SUCCESS) break; 

	}
      }


      // Create photons from the extended sources (clusters) and insert them
      // to the photon list.
      if (cluster_image!=NULL) {
	// Loop over all pixels of the the image:
	int xcount, ycount;
	double ra, dec;
	for(xcount=0; (xcount<cluster_image->width)&&(status==EXIT_SUCCESS); xcount++) {
	  for(ycount=0; (ycount<cluster_image->width)&&(status==EXIT_SUCCESS); ycount++) {
	    // Check whether the pixel lies CLOSE TO the FOV:
	    ra  = (xcount-cluster_image->width/2+0.5)*cluster_image->pixelwidth;
	    dec = (ycount-cluster_image->width/2+0.5)*cluster_image->pixelwidth;
	    struct vector v = unit_vector(ra, dec);

	    if (check_fov(&v, &telescope.nz, close_fov_min_align)==0) {
	      
	      // --- Generate Photons from the pixel.
	      
	      double random_number = get_random_number();
	      if (random_number < cluster_image->pixel[xcount][ycount].rate * dt *1.e9) {
		struct Photon new_photon; // buffer for new photon
		new_photon.ra  = ra;
		new_photon.dec = dec; 
		new_photon.direction = v; // REMOVE

		// Determine the energy of the new photon according to 
		// the default spectrum.
		new_photon.energy = photon_energy(spectrum_store.spectrum, detector);

		// Determine photon arrival time.
		double rnd_time = get_random_number();
		new_photon.time = time + dt*rnd_time;

		// Insert the photon into the time-ordered list.
		if ((status=insert_photon(&photon_list, new_photon))!=EXIT_SUCCESS) break;
	      }

	      // --- END of photon generation from the cluster image pixel.

	    } // END of check whether pixel is close to the FOV.
	  }
	} // END of loop over all pixel of the image.
      } // END  if(cluster_image!=NULL)



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
	if (check_fov(&photon_list->photon.direction,
		      &telescope.nz, fov_min_align)==0) {
	  // Photon is inside the FOV!

	  // Add the photon to the photon list file:
	  // Rescale from [rad] -> [deg]:
	  double ra  = photon_list->photon.ra *180./M_PI;
	  double dec = photon_list->photon.dec*180./M_PI;
	  fits_insert_rows(photonlist_fptr, photon_row++, 1, &status);
	  fits_write_col(photonlist_fptr, TDOUBLE, 1, photon_row, 1, 1, 
			 &photon_list->photon.time, &status);
	  fits_write_col(photonlist_fptr, TFLOAT, 2, photon_row, 1, 1, 
			 &photon_list->photon.energy, &status);
	  fits_write_col(photonlist_fptr, TDOUBLE, 3, photon_row, 1, 1, 
			 &ra, &status);
	  fits_write_col(photonlist_fptr, TDOUBLE, 4, photon_row, 1, 1, 
			 &dec, &status);
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
  
  headas_chat(5, "cleaning up ...\n");

  // Close FITS file
  if (photonlist_fptr) fits_close_file(photonlist_fptr, &status);

  // Release HEADAS random number generator:
  HDmtFree();
  gsl_rng_free(gsl_random_g);

  // Clear photon list
  clear_photon_list(&photon_list);

  // Release memory of orbit/attitude catalog
  if (sat_catalog) free(sat_catalog);

  // Point Sources
  free_PointSourceFiles(pointsourcefiles, &status);
  if (pointsourcecatalog != NULL) {
    if (pointsourcecatalog->sources != NULL) {
      int count;
      for(count=0; count<pointsourcecatalog->nsources; count++) {
	if (pointsourcecatalog->sources[count].lightcurve!=NULL) 
	  free(pointsourcecatalog->sources[count].lightcurve);
      }
      free(pointsourcecatalog->sources);
    }
    free (pointsourcecatalog);
  }

  // Cluster Images
  free_ClusterImage(cluster_image);
  
  // Release source spectra
  free_spectra(&spectrum_store, N_SPECTRA_FILES);

  // Detector data structure
  if (detector!=NULL) {
    if(detector->rmf!=NULL) free(detector->rmf);
    free(detector);
  }

  // String buffers for filenames of source files
  if (source_filename!=NULL) {
    int count;
    for(count=0; count<MAX_N_POINTSOURCEFILES; count++) {
      if(source_filename[count]!=NULL) free(source_filename[count]);
    }
    free(source_filename);
  }

  if (status==EXIT_SUCCESS) headas_chat(0, "finished successfully!\n\n");

  return(status);
}





