#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include "photon_generation.h"


/** Determines whether a given angle (in [rad]) lies within the
    specified range. The function returns a "1" if the angle lies
    within the specified range, otherwise the return value is "0". */
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



/** Insert photons from a photon list into a given photon list
    file. If the file already exists, the photons are inserted at the
    right positions in the file. Photons outside the FoV are
    neglected. */
int insertValidPhotonsIntoFile(PhotonListFile* plf,
			       struct PhotonOrderedListEntry** photon_list,  
			       AttitudeCatalog* ac,
			       struct Telescope telescope,
			       double t0, double timespan)
{
  int status=EXIT_SUCCESS;
  // Buffer for time-ordered photon list.
  struct PhotonOrderedListEntry* pl_entry;

  // Check if the photon list contains any photons.
  if (NULL==(*photon_list)) return(status);

  // Minimum cos-value for sources inside the FoV:
  // angle(telescope.nz,source) <= 1/2 * diameter
  const double fov_min_align = cos(telescope.fov_diameter/2.); 

  // Move the photon file pointer to the right entry for inserting
  // the photons from the list.
  double time=0.;
  int anynul=0;
  while (plf->row>0) {
    // Read the time from the file.
    fits_read_col(plf->fptr, TDOUBLE, plf->ctime, 
		  plf->row, 1, 1, &time, &time, &anynul, &status);
    if (status!=EXIT_SUCCESS) return(status);;

    if ((*photon_list)->photon.time >= time) break;
    
    plf->row--;
  }

  // SCAN PHOTON LIST to store the photons in a FITS file.
  while (NULL!=(*photon_list)) {
    // If all photons up to the actual time have been 
    // treated, break the loop.
    if ((*photon_list)->photon.time > t0+timespan) {
      break;
    }

    // Check whether the photon is inside the FoV:
    // First determine telescope pointing direction at the time of 
    // the photon arrival.
    telescope.nz = getTelescopePointing(ac, (*photon_list)->photon.time,
					&status);
    if (EXIT_SUCCESS!=status) break;

    // Compare the photon direction to the unit vector specifiing the 
    // direction of the telescope axis:
    if (check_fov(&(*photon_list)->photon.direction,
		  &telescope.nz, fov_min_align)==0) {
      // Photon is inside the FOV!

      // Add the photon to the photon list file:
      // Rescale from [rad] -> [deg]:
      double ra  = (*photon_list)->photon.ra *180./M_PI;
      double dec = (*photon_list)->photon.dec*180./M_PI;
      // Determine the right position to insert the photon:
      while (plf->row<plf->nrows) {
	// Read the time from the file.
	fits_read_col(plf->fptr, TDOUBLE, plf->ctime, 
		      plf->row+1, 1, 1, &time, &time, &anynul, &status);
	if (status!=EXIT_SUCCESS) return(status);;

	if ((*photon_list)->photon.time < time) break;
	
	plf->row++;
      }
      // Insert a new line with the photon data:
      fits_insert_rows(plf->fptr, plf->row++, 1, &status);
      fits_write_col(plf->fptr, TDOUBLE, plf->ctime, 
		     plf->row, 1, 1, &(*photon_list)->photon.time, &status);
      fits_write_col(plf->fptr, TFLOAT, plf->cenergy, 
		     plf->row, 1, 1, &(*photon_list)->photon.energy, &status);
      fits_write_col(plf->fptr, TDOUBLE, plf->cra, 
		     plf->row, 1, 1, &ra, &status);
      fits_write_col(plf->fptr, TDOUBLE, plf->cdec, 
		     plf->row, 1, 1, &dec, &status);
      plf->nrows++;
    } // END of: photon is inside the FOV?

    // Move to the next entry in the photon list and clear the current entry.
    pl_entry = (*photon_list)->next; // Buffer
    free(*photon_list);
    *photon_list = pl_entry;

    if (status!=EXIT_SUCCESS) break;

  } // END of scanning the photon list.

  return(status);
}



/* Generate photons for a catalog of point sources and a given
   telescope attitude. The photons lying in the FoV are inserted into
   the given photon list file using the function
   insertValidPhotonsIntoFile(). */
int createPhotonsFromPointSources(PhotonListFile* plf,
				  PointSourceFile* psf,
				  AttitudeCatalog* ac,
				  struct Telescope telescope,
				  struct RMF* rmf,
				  double t0, double timespan)
{
  int status = EXIT_SUCCESS;


  // Defines the mathematical meaning of 'close' in the context that for 
  // sources 'close to the FOV' the simulation creates a light curve.
  const double close_mult = 1.5; 
  // Minimum cos-value for point sources close to the FOV (in the direct
  // neighborhood).
  const double close_fov_min_align = cos(close_mult*telescope.fov_diameter/2.); 
  // Maximum cos-value (minimum sin-value) for sources within the
  // preselection band along the orbit. The width of the
  // preselection band is chosen to be twice the diameter of the
  // FoV. (angle(n,source) > 90-bandwidth)
  const double pre_max_align = sin(2*telescope.fov_diameter);


  // Normalized vector perpendicular to the orbital plane,
  // used for a preselection of point sources out of the
  // entire catalog.
  Vector preselection_vector;
  // Pointer source catalog.
  PointSourceCatalog* psc=NULL;
  // Normalized vector pointing in the direction of the source.
  Vector source_vector;
  // Counter for the different sources in the catalog.
  long source_counter;
  // Time-ordered photon list containing all created photons in the sky.
  struct PhotonOrderedListEntry* photon_list=NULL;  
  // Current time and time step.
  double time;
  const double dt = 0.1;
  // First run of the loop.
  int first = 1;

  // Loop over the given time interval.
  for (time=t0; time<t0+timespan; time+=dt) {
    headas_chat(0, "\rtime: %.3lf s ", time);
    fflush(NULL);

    // First determine telescope pointing direction at the current time.
    telescope.nz = getTelescopePointing(ac, time, &status);
    if (EXIT_SUCCESS!=status) break;


    // PRESELECTION of Point sources:
    // Preselection of sources from the comprehensive catalog to 
    // improve the performance of the simulation:
    if ((1==first)||
	(fabs(scalar_product(&preselection_vector, &telescope.nz)) > 
	 sin(0.5*telescope.fov_diameter))) {

      // Release old PointSourceCatalog:
      free_PointSourceCatalog(psc);
      psc = NULL;

      // Preselect sources from the entire source catalog according to the 
      // satellite's direction of motion.
      // Calculate normalized vector perpendicular to the orbit plane:
      preselection_vector = 
	normalize_vector(vector_product(telescope.nz,
					ac->entry[ac->current_entry].nx));
      psc=getPointSourceCatalog(psf, preselection_vector, pre_max_align, &status);
      if((EXIT_SUCCESS!=status)||(NULL==psc)) break;

    } // END of preselection

    // CREATE PHOTONS for all POINT sources CLOSE TO the FOV
    for (source_counter=0; source_counter<psc->nsources; source_counter++) {
      // Check whether the source is inside the FOV:
      // Compare the source direction to the unit vector specifiing the 
      // direction of the telescope:
      source_vector = 
	unit_vector(psc->sources[source_counter].ra, 
		    psc->sources[source_counter].dec);
      if (check_fov(&source_vector, &telescope.nz, close_fov_min_align) == 0) {
	    
	// The source is inside the FOV  => create photons:
	if ((status=create_PointSourcePhotons(&(psc->sources[source_counter]), 
					      time, dt, &photon_list, rmf))
	    !=EXIT_SUCCESS) break; 
	    
      } // END of check if source is close to the FoV
    } // End of loop over all sources in the catalog.
    if (EXIT_SUCCESS!=status) break;
    // END of photon generation.

    // Check the photon list and insert all photons inside the FoV
    // into the PhotonListFile.
    status=insertValidPhotonsIntoFile(plf,
				      &photon_list,  
				      ac,
				      telescope,
				      t0, timespan);
    if (EXIT_SUCCESS!=status) break;

    first = 0;
  } // END of the loop over the time interval.
  if (EXIT_SUCCESS!=status) return(status);

  // Clear photon list
  clear_PhotonList(&photon_list);

  return(status);
}



/* Generate photons for a catalog of point sources and a given
   telescope attitude. The photons lying in the FoV are inserted into
   the given photon list file using the function
   insertValidPhotonsIntoFile(). */
int createPhotonsFromExtendedSources(PhotonListFile* plf,
				     ExtendedSourceFile* esf,
				     AttitudeCatalog* ac,
				     struct Telescope telescope,
				     struct RMF* rmf,
				     double t0, double timespan)
{
  int status = EXIT_SUCCESS;

  // Defines the mathematical meaning of 'close' in the context that for 
  // sources 'close to the FOV' the simulation creates a light curve.
  const double close_mult = 1.5; 
  // Minimum cos-value for point sources close to the FOV (in the direct
  // neighborhood).
  const double close_fov_min_align = cos(close_mult*telescope.fov_diameter/2.); 
  // Maximum cos-value (minimum sin-value) for sources within the
  // preselection band along the orbit. The width of the
  // preselection band is chosen to be twice the diameter of the
  // FoV. (angle(n,source) > 90-bandwidth)
  const double pre_max_align = sin(2*telescope.fov_diameter);

  
  // Normalized vector perpendicular to the orbital plane,
  // used for a preselection of point sources out of the
  // entire catalog.
  Vector preselection_vector;
  // Pointer source catalog.
  ExtendedSourceCatalog* esc=NULL;
  // Normalized vector pointing in the direction of the source.
  Vector source_vector;
  // Counter for the different sources in the catalog.
  long source_counter;
  // Time-ordered photon list containing all created photons in the sky.
  struct PhotonOrderedListEntry* photon_list=NULL;  
  // Current time and time step.
  double time;
  const double dt = 0.1;
  // First run of the loop.
  int first = 1;

  // Loop over the given time interval.
  for (time=t0; time<t0+timespan; time+=dt) {
    headas_chat(0, "\rtime: %.3lf s ", time);
    fflush(NULL);

    // First determine telescope pointing direction at the current time.
    telescope.nz = getTelescopePointing(ac, time, &status);
    if (EXIT_SUCCESS!=status) break;


    // PRESELECTION of Point sources:
    // Preselection of sources from the comprehensive catalog to 
    // improve the performance of the simulation:
    if ((1==first)||
	(fabs(scalar_product(&preselection_vector, &telescope.nz)) > 
	 sin(0.5*telescope.fov_diameter))) {
      
      // Release old PointSourceCatalog:
      free_ExtendedSourceCatalog(esc);
      esc = NULL;

      // Preselect sources from the entire source catalog according to the 
      // satellite's direction of motion.
      // Calculate normalized vector perpendicular to the orbit plane:
      preselection_vector = 
	normalize_vector(vector_product(telescope.nz,
					ac->entry[ac->current_entry].nx));
      esc=getExtendedSourceCatalog(esf, preselection_vector, pre_max_align, &status);
      if((EXIT_SUCCESS!=status)||(NULL==esc)) break;
      
    } // END of preselection

    
    // CREATE PHOTONS for all EXTENDED sources CLOSE TO the FOV
    for (source_counter=0; source_counter<esc->nsources; source_counter++) {
      // Check whether the source is inside the FOV:
      // Compare the source direction to the unit vector specifying the 
      // direction of the telescope:
      source_vector = 
	unit_vector(esc->sources[source_counter].ra, 
		    esc->sources[source_counter].dec);
      if (check_fov(&source_vector, &telescope.nz, close_fov_min_align) == 0) {
	    
	// The source is inside the FOV  => create photons:
	if ((status=create_ExtendedSourcePhotons(&(esc->sources[source_counter]), 
						 time, dt, &photon_list, rmf))
	    !=EXIT_SUCCESS) break;
	    
      } // END of check if source is close to the FoV
    } // End of loop over all sources in the catalog.
    if (EXIT_SUCCESS!=status) break;
    // END of photon generation.

    // Check the photon list and insert all photons inside the FoV
    // into the PhotonListFile.
    status=insertValidPhotonsIntoFile(plf,
				      &photon_list,  
				      ac,
				      telescope,
				      t0, timespan);
    if (EXIT_SUCCESS!=status) break;
    
    first = 0;
  } // END of the loop over the time interval.
  if (EXIT_SUCCESS!=status) return(status);
  
  // Clear photon list
  clear_PhotonList(&photon_list);
  
  return(status);
}



/* Generate photons for a source image and a given telescope
   attitude. The photons lying in the FoV are inserted into the given
   photon list file using the function
   insertValidPhotonsIntoFile(). */
int createPhotonsFromSourceImage(PhotonListFile* plf,
				 SourceImage* si,
				 AttitudeCatalog* ac,
				 struct Telescope telescope,
				 struct RMF* rmf,
				 double t0, double timespan)
{
  int status = EXIT_SUCCESS;

  // Normalized vector pointing in the direction of the reference pixel.
  Vector refpixel_vector;
  // Normalized vector pointing in the direction of the source.
  Vector source_vector;
  // Counters for loops over all pixels of the SourceImage.
  int xcount, ycount;
  // Right ascension and declination of the individual pixels.
  double ra, dec; 
  // Buffer for currently chosen pixel.
  struct SourceImagePixel* pixel; 
  // Time-ordered photon binary tree.
  struct PhotonBinaryTreeEntry* photon_tree=NULL;
  // Time-ordered photon list containing all created photons in the sky.
  struct PhotonOrderedListEntry* photon_list=NULL;
  // Current time and time step.
  double time;
  const double dt = 1.0;
  // Buffer for random numbers.
  double rnd;
  // Buffer for new photons.
  Photon new_photon;

  // Defines the mathematical meaning of 'close' in the context that for 
  // sources 'close to the FOV' the simulation creates a light curve.
  const double close_mult = 1.5; 
  // Minimum cos-value for point sources close to the FOV (in the direct
  // neighborhood).
  const double close_fov_min_align = cos(close_mult*telescope.fov_diameter/2.); 


  // Loop over the given time interval.
  for (time=t0; time<t0+timespan; time+=dt) {
    headas_chat(0, "\rtime: %.3lf s ", time);
    fflush(NULL);

    // First determine telescope pointing direction at the current time.
    telescope.nz = getTelescopePointing(ac, time, &status);
    if (EXIT_SUCCESS!=status) break;

    // Create photons from the extended source image (clusters) and 
    // insert them to the photon list.

    // Vector in the direction of the reference pixel.
    refpixel_vector = unit_vector(si->crval1, si->crval2);

    // Check whether the the current telescope axis points to a direction 
    // close to the specified cluster image field (3°). (TODO: flexible value)
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
	// The reference pixel is neither at the north
	// nor at the south pole.
	k.x = -sin(si->crval1);
	k.y =  cos(si->crval1);
	    
	l.x = -sin(si->crval2) * cos(si->crval1);
	l.y = -sin(si->crval2) * sin(si->crval1);
	l.z =  cos(si->crval2);
      } // END reference pixel is at none of the poles.
	    
	    
      // Loop over all pixels of the the image:
      for(xcount=0; xcount<si->naxis1; xcount++) {
	for(ycount=0; ycount<si->naxis2; ycount++) {
		
	  // Check whether the pixel lies CLOSE TO the FOV:
		
	  // Vector in the direction of the current pixel.
	  source_vector.x = refpixel_vector.x + 
	    (xcount - si->crpix1 + 1.)*si->cdelt1 * k.x +
	    (ycount - si->crpix2 + 1.)*si->cdelt2 * l.x;
	  source_vector.y = refpixel_vector.y + 
	    (xcount - si->crpix1 + 1.)*si->cdelt1 * k.y +
	    (ycount - si->crpix2 + 1.)*si->cdelt2 * l.y;
	  source_vector.z = refpixel_vector.z + 
	    (xcount - si->crpix1 + 1.)*si->cdelt1 * k.z +
	    (ycount - si->crpix2 + 1.)*si->cdelt2 * l.z;
	  normalize_vector_fast(&source_vector);
		
	  if (check_fov(&source_vector, &telescope.nz, 
			close_fov_min_align)==0) {
	    
	    // Generate Photons for this pixel.

	    // Determine the right ascension and
	    // declination of the pixel.
	    calculate_ra_dec(source_vector, &ra, &dec);

	    // Current pixel.
	    pixel = &(si->pixel[xcount][ycount]);
		  
	    // TODO: needs to be replaced by exponential distribution:
	    rnd = sixt_get_random_number();
	    if (rnd < pixel->rate * dt) {
	      // Generate a new photon for this pixel
	      new_photon.ra=ra; 
	      new_photon.dec=dec;
	      new_photon.direction=source_vector;
	      new_photon.time = time + sixt_get_random_number()*dt;
	      
	      // Determine the energy of the new photon according to 
	      // the default spectrum.
	      new_photon.energy = 
		photon_energy(si->spectrumstore.spectrum, rmf);

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
	if (EXIT_SUCCESS!=status) break;
      } // END of loop over all pixel of the image.
      if (EXIT_SUCCESS!=status) break;
    } // END of check whether telescope axis points to the direction of the cluster field.

    // If a binary tree with photon entries is present, insert its entries to the 
    // time-ordered photon list.
    if (NULL!=photon_tree) {
      struct PhotonOrderedListEntry* photon_list_current = photon_list;
      if (EXIT_SUCCESS!=(status=CreateOrderedPhotonList(&photon_tree, &photon_list, 
							&photon_list_current))) break;
    }

    if (EXIT_SUCCESS!=status) break;
    // END of photon generation.


    // Check the photon list and insert all photons inside the FoV
    // into the PhotonListFile.
    status=insertValidPhotonsIntoFile(plf,
				      &photon_list,  
				      ac,
				      telescope,
				      t0, timespan);
    if (EXIT_SUCCESS!=status) break;

  } // END of time loop.
  if (EXIT_SUCCESS!=status) return(status);

  // Clear photon list
  clear_PhotonList(&photon_list);

  return(status);
}



/** Determine the first HDU in an X-ray source FITS file, determine
    the Source category, and set the file pointer to this HDU. */
int getFirstSourceHDU(fitsfile* fptr, SourceCategory* sc)
{
  int status=EXIT_SUCCESS;

  // Set to initial default value.
  *sc=INVALID_SOURCE;

  do { // Error handling loop 

    // Determine the number of HDUs in the FITS file and the current HDU.
    int n_hdus=0, hdu=0;
    if (fits_get_num_hdus(fptr, &n_hdus, &status)) break;
    fits_get_hdu_num(fptr, &hdu);
    // Determine the type of the first HDU (which should be IMAGE_HDU in 
    // any case).
    int hdu_type=0;
    if (fits_get_hdu_type(fptr, &hdu_type, &status)) break;
    headas_chat(5, " checking HDU %d/%d (HDU type %d) ...\n", 
		hdu, n_hdus, hdu_type);
      
    // Check whether the first (IMAGE) HDU is empty.
    // For this purpose check the header keyword NAXIS.
    int hdu_naxis=0;
    char comment[MAXMSG];
    if (fits_read_key(fptr, TINT, "NAXIS", &hdu_naxis, comment, 
		      &status)) break;
    headas_chat(5, " NAXIS: %d", hdu_naxis);
      
    if (0<hdu_naxis) {
      // The first (IMAGE) extension is not empty. So the source category 
      // is SOURCE_IMAGE.
      // Load the source image data from the HDU.
      headas_chat(5, " --> source type: SOURCE_IMAGES\n load data from current"
		  " HDU ...\n");
      *sc = SOURCE_IMAGES;
      
    } else {
      // The primary extension is empty, so move to the next extension.
      headas_chat(5, " --> empty\n move to next HDU ...\n");
      while (hdu<n_hdus) {
	if (fits_movrel_hdu(fptr, 1, &hdu_type, &status)) break;
	hdu++;
	headas_chat(5, " checking HDU %d/%d (HDU type %d) ...\n", 
		    hdu, n_hdus, hdu_type);

	if (IMAGE_HDU==hdu_type) {
	  headas_chat(5, " --> source type: SOURCE_IMAGES\n load data from "
		      "current HDU ...\n");
	  *sc = SOURCE_IMAGES;

	} else { 
	  // Assume BINARY_TBL (=2).

	  // Check whether it is a point or extended source catalog.
	  // Header keyword SRCTYPE must be either 'POINT' or 'EXTENDED'.
	  char comment[MAXMSG]; // String buffer.
	  char srctype[MAXMSG]; // String buffer for the source type.
	  if (fits_read_key(fptr, TSTRING, "SRCTYPE", srctype, 
			    comment, &status)) break;    
	  // Convert string to upper case. 
	  strtoupper(srctype);
	  // Check for the different possible source types.
	  if (0==strcmp(srctype, "POINT")) {
	    headas_chat(5, " --> source type: POINT_SOURCES\n load data from "
			"current HDU ...\n");
	    *sc = POINT_SOURCES;

	  } else if (0==strcmp(srctype, "EXTENDED")) {
	    headas_chat(5, " --> source type: EXTENDED_SOURCES\n load data from "
			"current HDU ...\n");
	    *sc = EXTENDED_SOURCES;

	  } else {
	    *sc = INVALID_SOURCE;
	  }

	}
      } // END of loop over the remaining HDUs after the primary extension.
      if (EXIT_SUCCESS!=status) break;

    } // END of check whether primary extension is empty.

  } while(0); // END of error handling loop.

  if (INVALID_SOURCE==*sc) {
    status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not find valid source data in input file!\n",
		   status);
    return(status);
  }

  return(status);
}



int photon_generation_main() 
{
  // Program parameters.
  struct Parameters parameters;
  
  // FITS file containing the input sources.
  fitsfile* sources_fptr=NULL;
  // Source category.
  SourceCategory sourceCategory=INVALID_SOURCE;
  // Data structure for point sources:
  PointSourceFile* psf=NULL;
  // Data structure for extended sources:
  ExtendedSourceFile* esf=NULL;
  // X-ray Cluster image:
  SourceImage* si=NULL;
  // WCS keywords from the input source file:
  struct wcsprm *wcs=NULL;
  // Number of WCS keyword sets.
  int nwcs=0; 

  // Catalog with attitude data.
  AttitudeCatalog* attitudecatalog=NULL;
  // Telescope data (like FOV diameter or focal length)
  struct Telescope telescope; 
  // RMF & EBOUNDS
  struct RMF* rmf=NULL;

  // Data structure for the output to the photon list FITS file.
  PhotonListFile photonlistfile={.fptr=NULL};

  int status = EXIT_SUCCESS;  // error status flag


  // Register HEATOOL
  set_toolname("photon_generation");
  set_toolversion("0.01");


  do {  // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    if((status=photon_generation_getpar(&parameters))) break;
    
    telescope.fov_diameter = parameters.fov_diameter; // in [rad]

    // Initialize HEADAS random number generator.
    HDmtInit(SIXT_HD_RANDOM_SEED);
    
    // Get the satellite catalog with the orbit and (telescope) attitude data:
    if (NULL==(attitudecatalog=get_AttitudeCatalog(parameters.attitude_filename,
						   parameters.t0, parameters.timespan, 
						   &status))) break;

    // Read RMF and EBOUNDS from the given file.
    rmf = loadRMF(parameters.rmf_filename, &status);
    if (EXIT_SUCCESS!=status) break;


    // Open the source file to check the contents.
    headas_chat(5, "open FITS file '%s' searching for X-ray sources ...\n",
		parameters.sources_filename);
    if (fits_open_file(&sources_fptr, parameters.sources_filename, READONLY, 
		       &status)) break;
    // Get the source type.
    status=getFirstSourceHDU(sources_fptr, &sourceCategory);
    if(EXIT_SUCCESS!=status) break;

    // Load the X-ray source data using the appropriate routines for the 
    // respective source category.
    if (POINT_SOURCES==sourceCategory) { 
      // Load the point source catalog from the current HDU.
      psf=get_PointSourceFile_fromHDU(sources_fptr, &status);
      if (status != EXIT_SUCCESS) break;
      
      // Read the header from the Point Source FITS file to obtain
      // the WCS keywords.
      char* header; // Buffer string for the FITS header.
      int nkeyrec, nreject;
      if (fits_hdr2str(sources_fptr, 1, NULL, 0, &header, &nkeyrec, &status)) break;
      // Determine the WCS header keywords from the header string.
      if ((wcsbth(header, nkeyrec, 0, 3, WCSHDR_PIXLIST, NULL, &nreject, &nwcs, &wcs))) {
	status=EXIT_FAILURE;
	HD_ERROR_THROW("Error: could not read WCS keywords from FITS header!\n", status);
	break;
      }
      // Release memory allocated by fits_hdr2str().
      free(header);
      header=NULL;

    } else if (EXTENDED_SOURCES==sourceCategory) { 
      // Load the point source catalog from the current HDU.
      esf=get_ExtendedSourceFile_fromHDU(sources_fptr, &status);
      if (status != EXIT_SUCCESS) break;
      
      // Read the header from the Extended Source FITS file to obtain
      // the WCS keywords.
      char* header; // Buffer string for the FITS header.
      int nkeyrec, nreject;
      if (fits_hdr2str(sources_fptr, 1, NULL, 0, &header, &nkeyrec, &status)) break;
      // Determine the WCS header keywords from the header string.
      if ((wcsbth(header, nkeyrec, 0, 3, WCSHDR_PIXLIST, NULL, &nreject, &nwcs, &wcs))) {
	status=EXIT_FAILURE;
	HD_ERROR_THROW("Error: could not read WCS keywords from FITS header!\n", status);
	break;
      }
      // Release memory allocated by fits_hdr2str().
      free(header);
      header=NULL;

    } else if (SOURCE_IMAGES==sourceCategory) {

      // Read the cluster images from the current FITS HDU.
      si=get_SourceImage_fromHDU(sources_fptr, &status);
      if(EXIT_SUCCESS!=status) break;

      // Read the header from the LAST source image FITS file to obtain
      // the WCS keywords.
      char* header; // Buffer string for the FITS header.
      int nkeyrec, nreject;
      if (fits_hdr2str(sources_fptr, 1, NULL, 0, &header, &nkeyrec, &status)) break;
      // Determine the WCS header keywords from the header string.
      if ((wcspih(header, nkeyrec, 0, 3, &nreject, &nwcs, &wcs))) {
	status=EXIT_FAILURE;
	HD_ERROR_THROW("Error: could not read WCS keywords from FITS header!\n", status);
	break;
      }
      // Release memory allocated by fits_hdr2str().
      free(header);
      header=NULL;
      
    } // END of different source categories.

    if (0==nwcs) {
      headas_chat(1, "### Warning: source file contains no appropriate WCS header keywords!\n");
    }

    // Check, whether a new photon list file should be generated,
    // whether an existing file should be opened for inserting 
    // the newly generated photons.
    if (1==parameters.overwrite_photonlist) {
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
	// If WCS keywords have been read from the source catalog / images file,
	// store them also in the photon list.
	if (fits_update_key(photonlistfile.fptr, TDOUBLE, "REFXCRVL", &wcs->crval[0],
			    "", &status)) break;
	if (fits_update_key(photonlistfile.fptr, TDOUBLE, "REFYCRVL", &wcs->crval[1],
			    "", &status)) break;
      } else {
	// Otherwise, if no WCS keywords are specified so far, use the RA=0 and 
	// DEC=0 as default values.
	double wcsbuffer=0.;
	if (fits_update_key(photonlistfile.fptr, TDOUBLE, "REFXCRVL", &wcsbuffer,
			    "", &status)) break;
	if (fits_update_key(photonlistfile.fptr, TDOUBLE, "REFYCRVL", &wcsbuffer,
			    "", &status)) break;
      }
    } else {
      // Open an existing photon list file for inserting the newly
      // generated photons.
      status = openPhotonListFile(&photonlistfile, parameters.photonlist_filename, 
				  READWRITE);
      if (EXIT_SUCCESS!=status) break;
    }
    // --- End of Initialization ---



    // --- Beginning of Photon generation process ---
    // Start the actual photon generation (after loading required data):
    headas_chat(5, "start photon generation process ...\n");

    if (POINT_SOURCES==sourceCategory) {

      status=createPhotonsFromPointSources(&photonlistfile,
					   psf,
					   attitudecatalog,
					   telescope,
					   rmf,
					   parameters.t0, 
					   parameters.timespan);
      if (EXIT_SUCCESS!=status) break;

    } else if (EXTENDED_SOURCES==sourceCategory) {

      status=createPhotonsFromExtendedSources(&photonlistfile,
					      esf,
					      attitudecatalog,
					      telescope,
					      rmf,
					      parameters.t0, 
					      parameters.timespan);
      if (EXIT_SUCCESS!=status) break;

    } else if (SOURCE_IMAGES==sourceCategory) {
      status=createPhotonsFromSourceImage(&photonlistfile,
					  si,
					  attitudecatalog,
					  telescope,
					  rmf,
					  parameters.t0, 
					  parameters.timespan);
      if (EXIT_SUCCESS!=status) break;

    }
    // --- End of photon generation ---

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(3, "\ncleaning up ...\n");

  status += closePhotonListFile(&photonlistfile);

  // Release HEADAS random number generator:
  HDmtFree();

  // Attitude Catalog
  free_AttitudeCatalog(attitudecatalog);

  // Close the Source file.
  if (NULL!=sources_fptr) fits_close_file(sources_fptr, &status);

  // PointSourceFile.
  free_PointSourceFile(psf);

  // Extended Source Images.
  free_SourceImage(si);
  si=NULL;
  
  if (status==EXIT_SUCCESS) headas_chat(0, "finished successfully!\n\n");
  return(status);
}



int photon_generation_getpar(struct Parameters* parameters)
{
  char msg[MAXMSG];           // Error message buffer.
  int status = EXIT_SUCCESS;  // Error status flag.

  // Get the filename of the Attitude file (FITS file):
  if ((status = PILGetFname("attitude_filename", parameters->attitude_filename))) {
    HD_ERROR_THROW("Error reading the filename of the attitude file!\n", status);
  }

  // Get the filename of the detector redistribution file (FITS file)
  else if ((status = PILGetFname("rmf_filename", parameters->rmf_filename))) {
    HD_ERROR_THROW("Error reading the filename of the detector" 
		   "redistribution matrix file (RMF)!\n", status);
  }

  // Get the start time of the photon generation
  else if ((status = PILGetReal("t0", &parameters->t0))) {
    HD_ERROR_THROW("Error reading the 't0' parameter!\n", status);
  }

  // Get the timespan for the photon generation
  else if ((status = PILGetReal("timespan", &parameters->timespan))) {
    HD_ERROR_THROW("Error reading the 'timespan' parameter!\n", status);
  }

  // Get the diameter of the FOV (in arcmin)
  else if ((status = PILGetReal("fov_diameter", &parameters->fov_diameter))) {
    HD_ERROR_THROW("Error reading the diameter of the FOV!\n", status);
  }

  // Determine the name of the file that contains the input sources (either
  // a point source catalog, source images, or a FITS grouping extension listing
  // several input files).
  else if ((status = PILGetFname("sources_filename", parameters->sources_filename))) {
    HD_ERROR_THROW("Error reading the filename of the input sources!\n", status);
  }

  // Get the filename of the Photon-List file (FITS file):
  else if ((status = PILGetFname("photonlist_filename", 
			    parameters->photonlist_filename))) {
    sprintf(msg, "Error reading the filename of the output file for "
	    "the photon list!\n");
    HD_ERROR_THROW(msg, status);
  }

  // Determine whether an existing photon list file should be
  // overwritten or not.
  else if ((status = PILGetInt("overwrite_photonlist", 
			       &parameters->overwrite_photonlist))) {
    HD_ERROR_THROW("Error reading the overwrite parameter!\n", status);
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
  if (EXIT_SUCCESS!=status) return(status);
  // Set the photon list template file:
  strcat(parameters->photonlist_template, "/photonlist.tpl");

  // Convert angles from [arc min] to [rad].
  parameters->fov_diameter *= M_PI/(180.*60.);

  return(status);
}


