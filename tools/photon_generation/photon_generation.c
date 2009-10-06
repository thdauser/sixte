#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include "photon_generation.h"


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
  if ((POINT_SOURCES != parameters->source_category) && 
      (SOURCE_IMAGES != parameters->source_category)) {
    status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: unknown source category file!\n", status);
  }
  if (EXIT_SUCCESS!=status) return(status);


  // Determine the name of the file that contains a list with the filenames of 
  // the source files.
  if ((status = PILGetFname("sourcelist_filename", 
			    parameters->sourcelist_filename))) {
    HD_ERROR_THROW("Error reading the filename of the list of source catalogs!\n", status);
  }

  // Get the filename of the Photon-List file (FITS file):
  else if ((status = PILGetFname("photonlist_filename", 
			    parameters->photonlist_filename))) {
    sprintf(msg, "Error reading the filename of the output file for "
	    "the photon list!\n");
    HD_ERROR_THROW(msg, status);
  }
  if (EXIT_SUCCESS!=status) return(status);

  // Get the name of the FITS template directory.
  // First try to read it from the environment variable.
  // If the variable does not exist, read it from the PIL.
  char* buffer;
  if (NULL!=(buffer=getenv("SIXT_FITS_TEMPLATES"))) {
    strcpy(parameters->photonlist_template, buffer);
  } else {
    if ((status = PILGetFname("fits_templates", parameters->photonlist_template))) {
      HD_ERROR_THROW("Error reading the path of the FITS templates!\n", status);
      
    }
  }
  // Set the photon list template file:
  strcat(parameters->photonlist_template, "/photonlist.tpl");

  // Convert angles from [arc min] to [rad].
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
  
  // Data structures for point sources:
  PointSourceFileCatalog* pointsourcefilecatalog=NULL;
  PointSourceCatalog* pointsourcecatalog=NULL;
  // X-ray Cluster image:
  SourceImageCatalog* sic = NULL;
  // WCS keywords from the input source file:
  struct wcsprm *wcs;
  // Number of WCS keyword sets.
  int nwcs; 

  // Catalog with attitude data.
  AttitudeCatalog* attitudecatalog=NULL;

  // Telescope data (like FOV diameter or focal length)
  struct Telescope telescope; 

  // RMF & EBOUNDS
  struct RMF* rmf=NULL;

  // Time-ordered photon binary tree.
  struct PhotonBinaryTreeEntry* photon_tree=NULL;
  // Time-ordered photon list containing all created photons in the sky
  struct PhotonOrderedListEntry* photon_list=NULL;  
  // Data structure for the output to the photon list FITS file.
  PhotonListFile photonlistfile;

  // Time step for sky scanning loop
  double dt = 0.1;

  char msg[MAXMSG];           // error message buffer
  int status = EXIT_SUCCESS;  // error status flag


  // Register HEATOOL
  set_toolname("photon_generation");
  set_toolversion("0.01");


  do {  // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    if((status=photon_generation_getpar(&parameters))) break;
    
    telescope.fov_diameter = parameters.fov_diameter; // in [rad]


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


    // Initialize HEADAS random number generator.
    HDmtInit(1);

    
    // Get the satellite catalog with the orbit and (telescope) attitude data:
    if (NULL==(attitudecatalog=get_AttitudeCatalog(parameters.attitude_filename,
						   parameters.t0, parameters.timespan, 
						   &status))) break;

    // Read RMF and EBOUNDS from the given file.
    rmf = loadRMF(parameters.rmf_filename, &status);
    if (EXIT_SUCCESS!=status) break;



    if (POINT_SOURCES==parameters.source_category) { 

      // Load the source catalogs from the files:
      pointsourcefilecatalog = get_PointSourceFileCatalog();
      if (NULL==pointsourcefilecatalog) {
	status = EXIT_FAILURE;
	sprintf(msg, "Error: allocation of PointSourceFileCatalog failed!\n");
	HD_ERROR_THROW(msg, status);
	break;
      }

      // Open the cluster list file:
      FILE* pointsourcelist_fptr = fopen(parameters.sourcelist_filename, "r");
      if (NULL==pointsourcelist_fptr) {
	sprintf(msg, "Error: could not open the file containing the list "
		"of point source catalogs!\n");
	HD_ERROR_THROW(msg, status);
	break;
      } else {
	// Determine the number of lines (= number of extended source files).
	char line[MAXMSG];
	pointsourcefilecatalog->nfiles = 0;
	while (fgets(line, MAXMSG, pointsourcelist_fptr)) {
	  pointsourcefilecatalog->nfiles++;
	}

	if (0==pointsourcefilecatalog->nfiles) {
	  sprintf(msg, "### Warning: Point Source List File '%s' contains no data!\n",
		  parameters.sourcelist_filename);
	}

	// Allocate memory for the point source files.
	pointsourcefilecatalog->files = 
	  (PointSourceFile**)malloc(pointsourcefilecatalog->nfiles*sizeof(PointSourceFile*));
	if (NULL==pointsourcefilecatalog->files) {
	  status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: not enough memory to allocate PointSourceFiles!\n", status);
	  break;
	}	

	// Load all point source catalogs specified in the point source list file.
	// Set the file pointer back to the beginning of the file:
	fseek(pointsourcelist_fptr, 0, SEEK_SET);
	int file_counter = 0;
	while (fscanf(pointsourcelist_fptr, "%s\n", line)>0) {
	  // Add the specified PointSourceFile the the PointSourceFileCatalog:
	  pointsourcefilecatalog->files[file_counter] = 
	    get_PointSourceFile_fromFile(line, &status);
	  if (status != EXIT_SUCCESS) break;
	}
	if (status != EXIT_SUCCESS) break;
      
	// Close the Point Source List File.
	fclose(pointsourcelist_fptr);

	// Read the header from the LAST Point Source FITS file to obtain
	// the WCS keywords.
	fitsfile* catalog_fptr=NULL;
	char* header; // Buffer string for the FITS header.
	int nkeyrec, nreject;
	if (fits_open_image(&catalog_fptr, line, READONLY, &status)) break;
	if (fits_hdr2str(catalog_fptr, 1, NULL, 0, &header, &nkeyrec, &status)) break;
	fits_close_file(catalog_fptr, &status);
	if (status!=EXIT_SUCCESS) break;	
	// Determine the WCS header keywords from the header string.
	if ((wcsbth(header, nkeyrec, 0, 3, WCSHDR_PIXLIST, NULL, &nreject, &nwcs, &wcs))) {
	  status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: could not read WCS keywords from FITS header!\n", status);
	  break;
	}
      }

      // Use a short time interval for the orbit update:
      dt = 0.01;
      
    } else if (SOURCE_IMAGES==parameters.source_category) {
      // Read the cluster images from the specified FITS files.

      // Get a new SourceImageCatalog object.
      sic = get_SourceImageCatalog();
      if (NULL==sic) {
	status = EXIT_FAILURE;
	HD_ERROR_THROW("Error: allocation of SourceImageCatalog failed!\n", status);
	break;
      }

      // Open the cluster list file:
      FILE* sourceimagelist_fptr = fopen(parameters.sourcelist_filename, "r");
      if (NULL==sourceimagelist_fptr) {
	HD_ERROR_THROW("Error: could not open the file containing the list "
		       "with the source images!\n", status);
	break;
      } else {

	// Determine the number of lines (= number of extended source files).
	char line[MAXMSG];
	while (fgets(line, MAXMSG, sourceimagelist_fptr)) {
	  sic->nimages++;
	}

	if (sic->nimages<=0) {
	  status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: invalid number of sources files with "
			 "extended sources!\n", status);
	  break;
	}
	
	sic->images = (SourceImage**)malloc(sic->nimages*sizeof(SourceImage*));
	if (NULL==sic->images) {
	  status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: memory allocation for ClusterImageCatalog failed!\n", status);
	  break;
	}

	// Load all cluster image files specified in the cluster list file.
	// Set the file pointer back to the beginning of the file:
	fseek(sourceimagelist_fptr, 0, SEEK_SET);
	int image_counter=0;
	while (fscanf(sourceimagelist_fptr, "%s\n", line)>0) {
	  // Load the specified galaxy cluster image:
	  sic->images[image_counter++] = get_SourceImage_fromFile(line, &status);
	  if (status != EXIT_SUCCESS) break;
	} // END of loop over all file entries in the cluster list file
	if (status != EXIT_SUCCESS) break;

	// Close the cluster list file:
	fclose(sourceimagelist_fptr);

	// Read the header from the LAST source image FITS file to obtain
	// the WCS keywords.
	fitsfile* image_fptr=NULL;
	char* header; // Buffer string for the FITS header.
	int nkeyrec, nreject;
	if (fits_open_image(&image_fptr, line, READONLY, &status)) break;
	if (fits_hdr2str(image_fptr, 1, NULL, 0, &header, &nkeyrec, &status)) break;
	fits_close_file(image_fptr, &status);
	if (status != EXIT_SUCCESS) break;	
	// Determine the WCS header keywords from the header string.
	if ((wcspih(header, nkeyrec, 0, 3, &nreject, &nwcs, &wcs))) {
	  status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: could not read WCS keywords from FITS header!\n", status);
	  break;
	}
      }
      
      // Use a relatively large \Delta t for the time loop:
      dt = 1.0;

    } else {
      status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: wrong source category!\n", status);
    } // END of different source categories.
    if (EXIT_SUCCESS!=status) break;

    if (0==nwcs) {
      headas_chat(1, "### Warning: source file contains no appropriate WCS header keywords!\n");
    }

    // Generate new photon list FITS file for output of generated photons.
    status = openNewPhotonListFile(&photonlistfile, parameters.photonlist_filename, 
				   parameters.photonlist_template);
    if (EXIT_SUCCESS!=status) break;

    // Add important HEADER keywords to the photon list.
    if (fits_update_key(photonlistfile.fptr, TSTRING, "ATTITUDE", 
			parameters.attitude_filename,
			"name of the attitude FITS file", &status)) break;
    // WCS keywords.
    if (nwcs > 0) {
      if (fits_update_key(photonlistfile.fptr, TDOUBLE, "REFXCRVL", &wcs->crval[0],
			  "", &status)) break;
      if (fits_update_key(photonlistfile.fptr, TDOUBLE, "REFYCRVL", &wcs->crval[1],
			  "", &status)) break;
    }
    // --- End of Initialization ---



    // --- Beginning of Photon generation process ---

    // LOOP over all timesteps given the specified timespan from t0 to t0+timespan
    double time; // current time

    long attitude_counter=0; // counter for orbit readout loop
    long last_attitude_counter=0; // stores sat_counter of former repetition, 
                                  // so the searching loop
                                  // doesn't have to start at 0 every time.

    // normalized vector perpendicular to the orbital plane
    Vector preselection_vector;


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

      // Get the last attitude entry before 'time' (in order to interpolate 
      // the attitude at this time between the neighboring calculated values):
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
      telescope.nz = interpolate_vec(attitudecatalog->entry[attitude_counter].nz, 
				     attitudecatalog->entry[attitude_counter].time, 
				     attitudecatalog->entry[attitude_counter+1].nz, 
				     attitudecatalog->entry[attitude_counter+1].time, 
				     time);
      normalize_vector_fast(&telescope.nz);


      if (parameters.source_category==POINT_SOURCES) {

	// PRESELECTION of Point sources
	// Preselection of sources from the comprehensive catalog to 
	// improve the performance of the simulation:
	if ((1==first)||
	    (fabs(scalar_product(&preselection_vector, &telescope.nz)) > 
	     sin(0.5*telescope.fov_diameter))) {

	  // Release old PointSourceCatalog:
	  free_PointSourceCatalog(pointsourcecatalog);
	  pointsourcecatalog = NULL;

	  // Preselect sources from the entire source catalog according to the 
	  // satellite's direction of motion.
	  // Calculate normalized vector perpendicular to the orbit plane:
	  preselection_vector = normalize_vector(vector_product(normalize_vector(
	      attitudecatalog->entry[attitude_counter].nz),
	      normalize_vector(attitudecatalog->entry[attitude_counter].nx)));

	  pointsourcecatalog=get_PointSourceCatalog(pointsourcefilecatalog, 
						    preselection_vector, 
						    pre_max_align,
						    &status);
	  if((EXIT_SUCCESS!=status)||(NULL==pointsourcecatalog)) break;
	}
	// END of preselection


	// CREATE PHOTONS for all POINT sources CLOSE TO the FOV
	long source_counter;
	for (source_counter=0; source_counter<pointsourcecatalog->nsources; 
	     source_counter++) {
	  // Check whether the source is inside the FOV:
	  	  
	  // Compare the source direction to the unit vector specifiing the 
	  // direction of the telescope:
	  Vector source_vector = 
	    unit_vector(pointsourcecatalog->sources[source_counter].ra, 
			pointsourcecatalog->sources[source_counter].dec);
	  if (check_fov(&source_vector, &telescope.nz, close_fov_min_align) == 0) {
	    
	    // The source is inside the FOV  => create photons:
	    if ((status=create_photons(&(pointsourcecatalog->sources[source_counter]), 
				       time, dt, &photon_list, rmf))
		!=EXIT_SUCCESS) break; 
	    
	  }
	}

      } else if (SOURCE_IMAGES==parameters.source_category) {

	// Create photons from the extended sources (clusters) and insert them
	// to the photon list.
	int image_counter;
	for (image_counter=0; image_counter<sic->nimages; image_counter++) {

	  // Vector in the direction of the reference pixel.
	  Vector refpixel_vector = 
	    unit_vector(sic->images[image_counter]->crval1, 
			sic->images[image_counter]->crval2);

	  // Check whether the the current telescope axis points to a direction 
	  // close to the specified cluster image field (3°).
	  if (check_fov(&refpixel_vector, &telescope.nz, cos(3.*M_PI/180.) )==0) {
	    // Vector in the direction of the 1st image coordinate (right ascension).
	    Vector k = {0., 0., 0.};
	    // Vector in the direction of the 2nd image coordinate (declination).
	    Vector l = {0., 0., 0.};

	    // Determine a local coordinate system for the cluster image.
	    if (fabs(refpixel_vector.z-1.) < 1.e-6) {
	      // The reference pixel lies at the north pole.
	      k.y = 1.;
	      l.x = -1.;
	    } else if (fabs(refpixel_vector.z+1) < 1.e-6) {
	      // The reference pixel lies at the south pole.
	      k.y = 1.;
	      l.x = 1.;
	    } else {
	      // The reference pixel is neither at the north nor at the south pole.
	      k.x = -sin(sic->images[image_counter]->crval1);
	      k.y =  cos(sic->images[image_counter]->crval1);
	    
	      l.x = -sin(sic->images[image_counter]->crval2) * 
		cos(sic->images[image_counter]->crval1);
	      l.y = -sin(sic->images[image_counter]->crval2) *
		sin(sic->images[image_counter]->crval1);
	      l.z =  cos(sic->images[image_counter]->crval2);
	    } // END reference pixel is at none of the poles.
	    
	    
	    // Loop over all pixels of the the image:
	    int xcount, ycount;
	    for(xcount=0; 
		(xcount<sic->images[image_counter]->naxis1)&&(EXIT_SUCCESS==status);
		xcount++) {
	      for(ycount=0; 
		  (ycount<sic->images[image_counter]->naxis2)&&(EXIT_SUCCESS==status); 
		  ycount++) {
		
		// Check whether the pixel lies CLOSE TO the FOV:
		
		// Vector in the direction of the current pixel.
		Vector pixel_vector; 
		pixel_vector.x = refpixel_vector.x + 
		  (xcount - sic->images[image_counter]->crpix1 + 1.)
		  *sic->images[image_counter]->cdelt1 * k.x +
		  (ycount - sic->images[image_counter]->crpix2 + 1.)
		  *sic->images[image_counter]->cdelt2 * l.x;
		pixel_vector.y = refpixel_vector.y + 
		  (xcount - sic->images[image_counter]->crpix1 + 1.)
		  *sic->images[image_counter]->cdelt1 * k.y +
		  (ycount - sic->images[image_counter]->crpix2 + 1.)
		  *sic->images[image_counter]->cdelt2 * l.y;
		pixel_vector.z = refpixel_vector.z + 
		  (xcount - sic->images[image_counter]->crpix1 + 1.)
		  *sic->images[image_counter]->cdelt1 * k.z +
		  (ycount - sic->images[image_counter]->crpix2 + 1.)
		  *sic->images[image_counter]->cdelt2 * l.z;
		normalize_vector_fast(&pixel_vector);
		
		if (check_fov(&pixel_vector, &telescope.nz, close_fov_min_align)==0) {
		  
		  // --- Generate Photons from the pixel.

		  // Determine the right ascension and declination of the pixel.
		  double ra, dec;
		  calculate_ra_dec(pixel_vector, &ra, &dec);

		  // Current pixel.
		  struct SourceImagePixel* pixel = 
		    &(sic->images[image_counter]->pixel[xcount][ycount]);

		  
		  // TODO: needs to be replaced by exponential distribution:
		  double rnd = get_random_number();
		  if (rnd < pixel->rate * dt) {
		    // Generate a new photon for this pixel
		    Photon new_photon  = {
		      .ra=ra, .dec=dec, .direction=pixel_vector,
		      .time = time + get_random_number()*dt };
		    
		    // Determine the energy of the new photon according to 
		    // the default spectrum.
		    new_photon.energy = 
		      photon_energy(sic->images[image_counter]->spectrumstore.spectrum, rmf);

		    // Add the photon to the binary tree.
		    if ((status=insert_Photon2BinaryTree(&photon_tree, &new_photon))
			!=EXIT_SUCCESS) break;
		  }
		  /*
		  if (pixel->t_next_photon < time) {
		    pixel->t_next_photon = time;
		  }
		  
		  while (pixel->t_next_photon <= time + dt) {

		    // Determine photon arrival time.
		    pixel->t_next_photon += rndexp(1./(double)pixel->rate);

		    Photon new_photon = { // Buffer for the new photon.
		      .ra=ra, .dec=dec, .direction=pixel_vector,
		      .time = pixel->t_next_photon };

		    // Determine the energy of the new photon according to 
		    // the default spectrum.
		    new_photon.energy = photon_energy(spectrum_store.spectrum, rmf);
		    
		    // Add the photon to the binary tree.
		    if ((status=insert_Photon2BinaryTree(&photon_tree, &new_photon))
			!=EXIT_SUCCESS) break;

		  }		  
		  */
		  // END of photon generation from the SourceImagePixel.
		
		} else { // END of check whether pixel is close to the FOV.
		  // The source image pixel is far away from the telescope axis.
		  // Therefore we don't have to consider the directly neighboring 
		  // pixels and do a bigger step than to the nearest pixel.
		  ycount+=10;
		}
	      }
	    } // END of loop over all pixel of the image.
	  } // END of check whether telescope axis points to the direction of the cluster field.
	} // END of loop over different extend source images in the catalog.
      } // END of decision which source category.



      // If a binary tree with photon entries is present, insert its entries to the 
      // time-ordered photon list.
      if (NULL!=photon_tree) {
	struct PhotonOrderedListEntry* photon_list_current = photon_list;
	if (EXIT_SUCCESS!=(status=CreateOrderedPhotonList(&photon_tree, &photon_list, 
							  &photon_list_current))) break;
      }


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
	telescope.nz = interpolate_vec(attitudecatalog->entry[attitude_counter].nz, 
				       attitudecatalog->entry[attitude_counter].time, 
				       attitudecatalog->entry[attitude_counter+1].nz, 
				       attitudecatalog->entry[attitude_counter+1].time, 
				       photon_list->photon.time);
	normalize_vector_fast(&telescope.nz);


	// Compare the photon direction to the unit vector specifiing the 
	// direction of the telescope axis:
	if (check_fov(&photon_list->photon.direction,
		      &telescope.nz, fov_min_align)==0) {
	  // Photon is inside the FOV!

	  // Add the photon to the photon list file:
	  // Rescale from [rad] -> [deg]:
	  double ra  = photon_list->photon.ra *180./M_PI; // TODO: photon.ra should already be rad
	  double dec = photon_list->photon.dec*180./M_PI;
	  fits_insert_rows(photonlistfile.fptr, photonlistfile.row++, 1, &status);
	  fits_write_col(photonlistfile.fptr, TDOUBLE, photonlistfile.ctime, 
			 photonlistfile.row, 1, 1, &photon_list->photon.time, &status);
	  fits_write_col(photonlistfile.fptr, TFLOAT, photonlistfile.cenergy, 
			 photonlistfile.row, 1, 1, &photon_list->photon.energy, &status);
	  fits_write_col(photonlistfile.fptr, TDOUBLE, photonlistfile.cra, 
			 photonlistfile.row, 1, 1, &ra, &status);
	  fits_write_col(photonlistfile.fptr, TDOUBLE, photonlistfile.cdec, 
			 photonlistfile.row, 1, 1, &dec, &status);
	  photonlistfile.nrows++;
	} // END of: photon is inside the FOV?

	// Move to the next entry in the photon list and clear the current entry.
	struct PhotonOrderedListEntry* pl_entry = photon_list->next; // Buffer
	free(photon_list);
	photon_list = pl_entry;

      } // END of scanning the photon list.

      last_attitude_counter = attitude_counter;
      first = 0;
    } // END of outer time loop.

    // --- End of photon creation process ---

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(5, "\ncleaning up ...\n");

  status += closePhotonListFile(&photonlistfile);

  // Release HEADAS random number generator:
  HDmtFree();

  // Clear photon list
  clear_PhotonList(&photon_list);

  // Attitude Catalog
  free_AttitudeCatalog(attitudecatalog);

  // Point Sources
  free_PointSourceFileCatalog(pointsourcefilecatalog);
  free_PointSourceCatalog(pointsourcecatalog);

  // Extended Source Images
  free_SourceImageCatalog(sic);
  
  if (status==EXIT_SUCCESS) headas_chat(0, "finished successfully!\n\n");
  return(status);
}



