#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include "photon_generation.h"
#include "clusters.c"


struct Parameters {
  char attitude_filename[FILENAME_LENGTH];
  char spectrum_filename[N_SPECTRA_FILES][FILENAME_LENGTH];
  char rmf_filename[FILENAME_LENGTH];
  char sourceimagelist_filename[FILENAME_LENGTH];
  char photonlist_filename[FILENAME_LENGTH];

  SourceCategory source_category;

  double t0, timespan;
  double fov_diameter;
};


////////////////////////////
int photon_generation_getpar(struct Parameters* parameters)
{
  char msg[MAXMSG];           // error message buffer
  int status = EXIT_SUCCESS;  // error status flag


  // Get the filename of the Attitude file (FITS file):
  if ((status = PILGetFname("attitude_filename", parameters->attitude_filename))) {
    sprintf(msg, "Error reading the filename of the attitude file!\n");
    HD_ERROR_THROW(msg, status);
  }

  // Get the filename of the default source spectrum (PHA FITS file)
  if ((status = PILGetFname("spectrum_filename", parameters->spectrum_filename[0]))) {
    sprintf(msg, "Error reading the filename of the default spectrum (PHA)!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the filename of the detector redistribution file (FITS file)
  else if ((status = PILGetFname("rmf_filename", parameters->rmf_filename))) {
    sprintf(msg, "Error reading the filename of the detector" 
	    "redistribution matrix file (RMF)!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the start time of the photon generation
  else if ((status = PILGetReal("t0", &parameters->t0))) {
    sprintf(msg, "Error reading the 't0' parameter!\n");
    HD_ERROR_THROW(msg, status);
  }

  // Get the timespan for the photon generation
  else if ((status = PILGetReal("timespan", &parameters->timespan))) {
    sprintf(msg, "Error reading the 'timespan' parameter!\n");
    HD_ERROR_THROW(msg, status);
  }

  // Get the diameter of the FOV (in arcmin)
  else if ((status = PILGetReal("fov_diameter", &parameters->fov_diameter))) {
    sprintf(msg, "Error reading the diameter of the FOV!\n");
    HD_ERROR_THROW(msg, status);
  }

  // Determine the category of the input sources:
  int category;
  if ((status = PILGetInt("source_category", &category))) {
    sprintf(msg, "Error: wrong source category!\n");
    HD_ERROR_THROW(msg, status);
  }
  parameters->source_category=category;
  if (EXIT_SUCCESS!=status) return(status);

  if (EXTENDED_SOURCES==parameters->source_category) {
    // Determine the name of the file that contains the filenames of 
    // the extended source files.
    if ((status = PILGetFname("sourceimagelist_filename", 
			      parameters->sourceimagelist_filename))) {
      sprintf(msg, "Error reading the filename of the cluster list file!\n");
      HD_ERROR_THROW(msg, status);
    }
  }

  // Get the filename of the Photon-List file (FITS file):
  if ((EXIT_SUCCESS==status) && 
      (status = PILGetFname("photonlist_filename", 
			    parameters->photonlist_filename))) {
    sprintf(msg, "Error reading the filename of the output file for "
	    "the photon list!\n");
    HD_ERROR_THROW(msg, status);
  }

  // Convert angles from [arc min] to [rad]
  parameters->fov_diameter = parameters->fov_diameter*M_PI/(180.*60.);

  return(status);
}



//////////////////////////
/** Determines whether a given angle (in [rad]) lies within the specified range. 
 * The function returns a "1" if the angle lies within the specified range, 
 * otherwise the return value is "0". */
int check_angle_range(double angle, double min, double max) 
{
  while (fabs(angle-min) > M_PI) { min-=2*M_PI; }
  while (fabs(angle-max) > M_PI) { max+=2*M_PI; }

  if ((angle>min)&&(angle<max)) {
    return(1);
  } else {
    return(0);
  }
}




//////////////////////////
int photon_generation_main() 
{
  // Program parameters.
  struct Parameters parameters;
  
  // Several input source catalog files:
  int n_sourcefiles;   // number of input source files
  // Filenames of the individual source catalog files (FITS):
  char** source_filename=NULL;
  // X-ray Cluster image:
  ClusterImageCatalog* cic = NULL;

  // New data structures for point sources:
  PointSourceFiles* pointsourcefiles=NULL;
  PointSourceCatalog* pointsourcecatalog=NULL;
  long source_counter;

  // Catalog with orbit and attitude data over a particular timespan
  //  struct Telescope *sat_catalog=NULL;     
  // Number of entries in the orbit list ( <= orbit_nrows)
  //  long sat_nentries;                      
  AttitudeCatalog* attitudecatalog=NULL;

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

  // time step for sky scanning loop
  double dt = 0.1;

  gsl_rng *gsl_random_g=NULL; // pointer to GSL random number generator

  char msg[MAXMSG];           // error message buffer
  int status = EXIT_SUCCESS;  // error status flag


  // Register HEATOOL
  set_toolname("photon_generation");
  set_toolversion("0.01");


  do {  // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    detector = get_Detector(&status);
    if(status!=EXIT_SUCCESS) break;
    

    if((status=photon_generation_getpar(&parameters))) break;
    
    telescope.fov_diameter = parameters.fov_diameter;

    // Set last_update to such a small value, that a preselection of the 
    // source catalog is performed at the first timestep (last_update 
    // contains the time of the last source catalog preselection.):
    //    const double former_time = parameters.t0 - ORBIT_UPDATE_TIME - 100.;
    //    double last_update = former_time;

    // Defines the mathematical meaning of 'close' in the context that for 
    // sources 'close to the FOV' the simulation creates a light curve.
    const double close_mult = 1.5; 
  
    // Calculate the minimum cos-value for sources inside the FOV and 
    // in the preselection band respectively:
    
    // minimum cos-value for sources inside the FOV 
    // angle(x0,source) <= 1/2 * diameter
    const double fov_min_align = cos(telescope.fov_diameter/2.); 
    
    /** Minimum cos-value for sources close to the FOV (in the direct neighborhood). */
    const double close_fov_min_align = cos(close_mult*telescope.fov_diameter/2.); 

    /** Maximum cos-value (minimum sin-value) for sources within the 
     * preselection band along the orbit. The width of the preselection band
     * is chose to be twice the diameter of the FOV. 
     * (angle(n,source) > 90-bandwidth) */
    const double pre_max_align = sin(2*telescope.fov_diameter);


    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);
    gsl_rng_env_setup();
    gsl_random_g = gsl_rng_alloc(gsl_rng_default);

    
    // Get the satellite catalog with the orbit and (telescope) attitude data:
    if (NULL==(attitudecatalog=get_AttitudeCatalog(parameters.attitude_filename,
						   parameters.t0, parameters.timespan, 
						   &status))) break;


    // Read the detector RMF and EBOUNDS from the specified file and assign them to the 
    // Detector data structure.
    if ((status=detector_assign_rsp(detector, parameters.rmf_filename)) 
	!= EXIT_SUCCESS) break;

    // Get the source spectra:
    if ((status=get_spectra(&spectrum_store, detector->rmf->NumberChannels, 
			    parameters.spectrum_filename, N_SPECTRA_FILES)) 
	!= EXIT_SUCCESS) break;



    if (parameters.source_category==POINT_SOURCES) { 

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

      // Get the number of source input-files
      if ((status = PILGetInt("n_sourcefiles", &n_sourcefiles))) {
	sprintf(msg, "Error reading the number of source catalog files!\n");
	HD_ERROR_THROW(msg, status);
	break;
      }
      // Get the filenames of the individual source catalogs.
      else {
	int counter;
	// filename-buffer to access the different source files
	char cbuffer[FILENAME_LENGTH];

	for(counter=0; counter<MAX_N_POINTSOURCEFILES; counter++) {
	  sprintf(cbuffer,"sourcefile%d", counter+1);

	  if (counter<n_sourcefiles) {
	    // Read the source file from using PIL.
	    if ((status = PILGetFname(cbuffer, source_filename[counter]))) {
	      sprintf(msg, "Error reading the name of the sourcefile No %d!\n", 
		      counter);
	      HD_ERROR_THROW(msg, status);
	    }
	  } else {
	    // Fill redundant input slots for source files with NULL value,
	    // in order to avoid errors with the HD_PARSTAMP routine.
	    PILPutFname(cbuffer, "");
	  }
	}
	// Clear the PIL parameter for the cluster filename in order to 
	// avoid problems with HD_PARSTAMP.
	PILPutFname("cluster_filename", "");

	if (status) break;
      }

      // Load the source catalogs from the files:
      pointsourcefiles = get_PointSourceFiles(n_sourcefiles, source_filename, &status);
      if (status != EXIT_SUCCESS) break;

      // Use a short time interval for the orbit update:
      dt = 0.001;
      
    } else if (parameters.source_category==EXTENDED_SOURCES) {
      // Read the cluster images from the specified FITS files.
      cic = get_ClusterImageCatalog();

      // Open the cluster list file:
      FILE* sourceimagelist_fptr = fopen(parameters.sourceimagelist_filename, "r");
      if (NULL==sourceimagelist_fptr) {
	sprintf(msg, "Error: could not open the file containing the list "
		"with the source images!\n");
	HD_ERROR_THROW(msg, status);
	break;
      } else {
	// Determine the number of lines (= number of extended source files).
	char line[MAXMSG];
	while (fgets(line, MAXMSG, sourceimagelist_fptr)) {
	  cic->nimages++;
	}

	if (cic->nimages<=0) {
	  status=EXIT_SUCCESS;
	  sprintf(msg, "Error: invalid number of sources files with "
		  "extended sources!\n");
	  HD_ERROR_THROW(msg, status);
	  break;
	}
	
	cic->images = (ClusterImage**)malloc(cic->nimages*sizeof(ClusterImage));
	if (NULL==cic->images) {
	  status=EXIT_SUCCESS;
	  sprintf(msg, "Error: memory allocation for ClusterImageCatalog failed!\n");
	  HD_ERROR_THROW(msg, status);
	  break;
	}

	// Load all cluster image files specified in the cluster list file.
	// Set the file pointer back to the beginning of the file:
	fseek(sourceimagelist_fptr, 0, SEEK_SET);
	int image_counter=0;
	while (fscanf(sourceimagelist_fptr, "%s\n", line)>0) {
	  // Load the specified galaxy cluster image:
	  cic->images[image_counter++] = get_ClusterImage_fromFile(line, &status);
	  if (status != EXIT_SUCCESS) break;
	} // END of loop over all file entries in the cluster list file

	// Close the cluster list file:
	fclose(sourceimagelist_fptr);
      }
      
      // Clear the filenames of the point source catalogs in the PIL parameter file,
      // otherwise HD_PARSTAMP might cause an error.

      // Filename-buffer to access the different source files:
      int file_counter;
      char cbuffer[FILENAME_LENGTH];
      for(file_counter=0; file_counter<MAX_N_POINTSOURCEFILES; file_counter++) {
	sprintf(cbuffer,"sourcefile%d", file_counter+1);
	// Fill redundant input slots for source files with NULL value,
	// in order to avoid errors with the HD_PARSTAMP routine.
	PILPutFname(cbuffer, "");
      }

      // Use a relative large \Delta t for the time loop:
      dt = 1.0;

    } else {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: wrong source category!\n");
      HD_ERROR_THROW(msg, status);
    } // different source categories


    // Delete old photon list FITS file:
    remove(parameters.photonlist_filename);
    // Create a new FITS file for the output of the photon list:
    if ((create_photonlist_file(&photonlist_fptr, 
				parameters.photonlist_filename, &status)))
      break;

    // Add important HEADER keywords to the photon list
    if (fits_write_key(photonlist_fptr, TSTRING, "ATTITUDE", 
		       parameters.attitude_filename,
		       "name of the attitude FITS file", &status)) break;
    


    // --- End of Initialization ---



    // --- Beginning of Photon generation process ---

    // LOOP over all timesteps given the specified timespan from t0 to t0+timespan
    double time;             // current time

    long attitude_counter=0;      // counter for orbit readout loop
    long last_attitude_counter=0; // stores sat_counter of former repetition, 
                                  // so the searching loop
                                  // doesn't have to start at 0 every time.
    long photon_row=0;       // current row in photon list FITS file;

    // normalized vector perpendicular to the orbital plane
    struct vector preselection_vector;

    struct Photon_Entry *pl_entry=NULL; // "counter" variable for the photon list


    // Beginning of actual simulation (after loading required data):
    headas_chat(5, "start photon generation process ...\n");


    // Time loop:
    // Timesteps are typically a fraction (e.g. 1/10) of the time, the satellite 
    // takes to slew over the entire FOV.
    int first = 1;  // first run of the loop
    for(time=parameters.t0; 
	(time<parameters.t0+parameters.timespan)&&(status==EXIT_SUCCESS); time+=dt) {
      headas_chat(0, "\rtime: %.3lf s ", time);
      fflush(NULL);

      // Get the last attitude entry before 'time'
      // (in order to interpolate the attitude at this time between 
      // the neighboring calculated values):
      for( ; attitude_counter<attitudecatalog->nentries-1; attitude_counter++) {
	if(attitudecatalog->entry[attitude_counter+1].time>time) {
	  break;
	}
      }
      if(fabs(attitudecatalog->entry[attitude_counter].time-time)>600.) { 
	// no entry within 10 minutes !!
	status = EXIT_FAILURE;
	sprintf(msg, "Error: no adequate orbit entry for time %lf!\n", time);
	HD_ERROR_THROW(msg,status);
	break;
      }


      // First determine telescope pointing direction at the current time.
      // TODO: replace this calculation by proper attitude interpolation.
      telescope.nz = 
	normalize_vector(interpolate_vec(attitudecatalog->entry[attitude_counter].nz, 
					 attitudecatalog->entry[attitude_counter].time, 
					 attitudecatalog->entry[attitude_counter+1].nz, 
					 attitudecatalog->entry[attitude_counter+1].time, 
					 time));

      if (parameters.source_category==POINT_SOURCES) {

	/*
	// PRESELECTION of Point sources
	// Preselection of sources from the comprehensive catalog to 
	// improve the performance of the simulation:
	if (attitudecatalog->entry[attitude_counter].time-last_update > ORBIT_UPDATE_TIME) {
	  // Preselect sources from the entire source catalog according to the 
	  // satellite's direction of motion.
	  // Calculate normalized vector perpendicular to the orbit plane:
	  preselection_vector = normalize_vector(vector_product(normalize_vector(
	      attitudecatalog->entry[attitude_counter].nz), 
	      normalize_vector(attitudecatalog->entry[attitude_counter].nx)));
	  if((status=get_PointSourceCatalog(pointsourcefiles, &pointsourcecatalog,
					    preselection_vector, pre_max_align,
					    spectrum_store))
	     !=EXIT_SUCCESS) break;
	  
	  // Update the catalog-update-counter
	  last_update = attitudecatalog->entry[attitude_counter].time;
	}
	// END of preselection
	*/
	// PRESELECTION of Point sources
	// Preselection of sources from the comprehensive catalog to 
	// improve the performance of the simulation:
	if ((1==first)||
	    (scalar_product(&preselection_vector, &telescope.nz) > 
	     cos(M_PI/2 - parameters.fov_diameter))) {
	  // Preselect sources from the entire source catalog according to the 
	  // satellite's direction of motion.
	  // Calculate normalized vector perpendicular to the orbit plane:
	  preselection_vector = normalize_vector(vector_product(normalize_vector(
	      attitudecatalog->entry[attitude_counter].nz), 
	      normalize_vector(attitudecatalog->entry[attitude_counter].nx)));
	  if((status=get_PointSourceCatalog(pointsourcefiles, &pointsourcecatalog,
					    preselection_vector, pre_max_align,
					    spectrum_store))
	     !=EXIT_SUCCESS) break;
	  
	  // Update the catalog-update-counter
	  //	  last_update = attitudecatalog->entry[attitude_counter].time;
	}
	// END of preselection


	// CREATE PHOTONS for all POINT sources CLOSE TO the FOV
	for (source_counter=0; source_counter<pointsourcecatalog->nsources; 
	     source_counter++) {
	  // Check whether the source is inside the FOV:
	  	  
	  // Compare the source direction to the unit vector specifiing the 
	  // direction of the telescope:
	  struct vector source_vector = 
	    unit_vector(pointsourcecatalog->sources[source_counter].ra, 
			pointsourcecatalog->sources[source_counter].dec);
	  if (check_fov(&source_vector, &telescope.nz, close_fov_min_align) == 0) {
	    
	    // The source is inside the FOV  => create photons:
	    if ((status=create_photons(&(pointsourcecatalog->sources[source_counter]), 
				       time, dt, &photon_list, detector, gsl_random_g))
		!=EXIT_SUCCESS) break; 
	    
	  }
	}

      } else if (EXTENDED_SOURCES==parameters.source_category) {

	// Create photons from the extended sources (clusters) and insert them
	// to the photon list.
	int image_counter;
	for (image_counter=0; image_counter<cic->nimages; image_counter++) {

	  // TODO
	  // Check whether the the current telescope axis lies within the specified field
	  // or CLOSE TO it. !!

	  // Loop over all pixels of the the image:
	  int xcount, ycount;
	  double ra, dec;
	  for(xcount=0; 
	      (xcount<cic->images[image_counter]->naxis1)&&(status==EXIT_SUCCESS); 
	      xcount++) {
	    for(ycount=0; 
		(ycount<cic->images[image_counter]->naxis2)&&(status==EXIT_SUCCESS); 
		ycount++) {
	      // Check whether the pixel lies CLOSE TO the FOV:
	      ra=cic->images[image_counter]->crval1+
		(xcount-cic->images[image_counter]->crpix1+0.5)
		*cic->images[image_counter]->cdelt1; // [rad]
	      dec=cic->images[image_counter]->crval2+
		(ycount-cic->images[image_counter]->crpix2+0.5)
		*cic->images[image_counter]->cdelt2; // [rad]
	      struct vector v = unit_vector(ra, dec);

	      if (check_fov(&v, &telescope.nz, close_fov_min_align)==0) {
	      
		// --- Generate Photons from the pixel.
		
		double random_number = get_random_number();
		if(random_number <
		   cic->images[image_counter]->pixel[xcount][ycount].rate*dt*4){
		  struct Photon new_photon = { // buffer for new photon
		    .ra=ra, .dec=dec, .direction=v }; 
		  
		  // Determine the energy of the new photon according to 
		  // the default spectrum.
		  new_photon.energy = photon_energy(spectrum_store.spectrum, detector);
		  
		  // Determine photon arrival time.
		  double rnd_time = get_random_number();
		  new_photon.time = time + dt*rnd_time;
		  
		  // Insert the photon into the time-ordered list.
		  if ((status=insert_photon(&photon_list, new_photon))!=EXIT_SUCCESS) 
		    break;
		}		
		// --- END of photon generation from the cluster image pixel.
		
	      } else { // END of check whether pixel is close to the FOV.
		// The source image pixel is far away from the telescope axis.
		// Therefore we don't have to consider the directly neighboring 
		// pixels and do a bigger step than to the nearest pixel.
		ycount+=9;
	      }
	    }
	  } // END of loop over all pixel of the image.
	} // END of loop over different extend source images in the catalog
      } // END of decision which source category


      // SCAN PHOTON LIST
      // Perform the actual measurement (i.e., scan the photon list).
      attitude_counter = last_attitude_counter; 
      // Because we are now searching for events in the time interval [time-dt,time].
      while ((photon_list != NULL) && (status == EXIT_SUCCESS)) {
	// If all photons up to the actual time have been treated, break the loop.
	if ((photon_list->photon.time > time + dt)||
	    (photon_list->photon.time > parameters.t0+parameters.timespan)) {
	  break;
	}

	// Get the last orbit entry before the time 'photon_list->photon.time'
	// (in order to interpolate the position and velocity at this time  between 
	// the neighboring calculated orbit positions):
	for( ; attitude_counter<attitudecatalog->nentries-1; attitude_counter++) {
	  if(attitudecatalog->entry[attitude_counter+1].time>photon_list->photon.time) {
	    break;
	  }
	}
	if(fabs(attitudecatalog->entry[attitude_counter].time-time)>600.) { 
	  // no entry within 10 minutes !!
	  status = EXIT_FAILURE;
	  sprintf(msg, "Error: no adequate orbit entry for time %lf!\n", time);
	  HD_ERROR_THROW(msg,status);
	  break;
	}

	// Check whether the photon is inside the FOV:
	// First determine telescope pointing direction at the current time.
	telescope.nz = 
	  normalize_vector(interpolate_vec(attitudecatalog->entry[attitude_counter].nz, 
					   attitudecatalog->entry[attitude_counter].time, 
					   attitudecatalog->entry[attitude_counter+1].nz, 
					   attitudecatalog->entry[attitude_counter+1].time, 
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
	} // END of photon is inside the FOV

	// Move to the next entry in the photon list and clear the current entry.
	pl_entry = photon_list->next_entry;
	free(photon_list);
	photon_list = pl_entry;

      }  // END of scanning the photon list.

      last_attitude_counter = attitude_counter;
      first = 0;
    }   // END of outer time loop.

    // --- End of photon creation process ---

  } while(0); // END of ERROR HANDLING Loop.



  // --- Clean up ---
  
  headas_chat(5, "\ncleaning up ...\n");

  // Close FITS file
  if (photonlist_fptr) fits_close_file(photonlist_fptr, &status);

  // Release HEADAS random number generator:
  HDmtFree();
  if (gsl_random_g!=NULL) gsl_rng_free(gsl_random_g);

  // Clear photon list
  clear_photon_list(&photon_list);

  // Release memory of orbit/attitude catalog
  //  if (sat_catalog) free(sat_catalog);

  free_AttitudeCatalog(attitudecatalog);

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

  // ClusterImageCatalog
  free_ClusterImageCatalog(cic);
  
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





