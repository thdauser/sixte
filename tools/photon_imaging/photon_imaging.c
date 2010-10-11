#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include "photon_imaging.h"


////////////////////////////////////
/** Main procedure. */
int photon_imaging_main() {
  struct Parameters parameters;

  AttitudeCatalog* ac=NULL;

  PhotonListFile photonlistfile;
  ImpactListFile* impactlistfile=NULL;
  // WCS reference values for the position of the orginial input data.
  // These data are needed for the eROSITA image reconstruction algorithm
  // in order to determine the right WCS header keywords for, e.g., cluster images.
  double refxcrvl, refycrvl; 

  // Detector data structure including telescope information like the PSF,
  // vignetting function, focal length, and FOV diameter.
  GenDet* det=NULL;

  struct Telescope telescope;

  int status=EXIT_SUCCESS; // Error status


  // Register HEATOOL:
  set_toolname("photon_imaging");
  set_toolversion("0.01");


  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // read parameters using PIL library
    if ((status=photon_imaging_getpar(&parameters))) break;

    // Initialize the detector data structure.
    det = newGenDet(parameters.xml_filename, &status);
    if (EXIT_SUCCESS!=status) break;

    // Calculate the minimum cos-value for sources inside the FOV: 
    // (angle(x0,source) <= 1/2 * diameter)
    const double fov_min_align = cos(det->fov_diameter/2.); 
    
    // Initialize HEADAS random number generator. Add the telescope number
    // to the standard seed in order to avoid getting the same random number
    // sequence for all the 7 telescopes. This results in problems due to
    // event correlation.
    HDmtInit(SIXT_HD_RANDOM_SEED+parameters.telescope);

    // Open the FITS file with the input photon list:
    status = openPhotonListFile(&photonlistfile, parameters.photonlist_filename, 
				READONLY);
    if (EXIT_SUCCESS!=status) break;
    // Read WCS keywords from FITS header.
    char comment[MAXMSG]; // String buffer.
    if (fits_read_key(photonlistfile.fptr, TDOUBLE, "REFXCRVL", &refxcrvl, 
		      comment, &status)) break;    
    if (fits_read_key(photonlistfile.fptr, TDOUBLE, "REFYCRVL", &refycrvl, 
		      comment, &status)) break;    

    // Open the attitude file specified in the header keywords of the photon list.
    if (fits_read_key(photonlistfile.fptr, TSTRING, "ATTITUDE", 
		      &parameters.attitude_filename, 
		      comment, &status)) break;
    if (NULL==(ac=get_AttitudeCatalog(parameters.attitude_filename,
				      parameters.t0, parameters.timespan, 
				      &status))) break;

    // Create a new FITS file for the output of the impact list.
    impactlistfile = openNewImpactListFile(parameters.impactlist_filename, 
					   parameters.impactlist_template,
					   &status);
    if (EXIT_SUCCESS!=status) break;

    // Write WCS header keywords.
    if (fits_update_key(impactlistfile->fptr, TDOUBLE, "REFXCRVL", &refxcrvl, "", &status)) break;
    if (fits_update_key(impactlistfile->fptr, TDOUBLE, "REFYCRVL", &refycrvl, "", &status)) break;
    // Add attitude filename.
    if (fits_update_key(impactlistfile->fptr, TSTRING, "ATTITUDE", parameters.attitude_filename,
		       "name of the attitude FITS file", &status)) break;
    
    // --- END of Initialization ---



    // --- Beginning of Imaging Process ---

    // Beginning of actual simulation (after loading required data):
    headas_chat(5, "start imaging process ...\n");

    // SCAN PHOTON LIST    
    // LOOP over all timesteps given the specified timespan from t0 to t0+timespan
    Photon photon={.time=0.};
    // Buffer for impact position:
    struct Point2d position;
    // Buffer for scalar product:
    double scp;
    while (photonlistfile.row<photonlistfile.nrows) {

      if (EXIT_SUCCESS!=status) break;
      
      // Read an entry from the photon list:
      status = PhotonListFile_getNextRow(&photonlistfile, &photon);
      if (status!=EXIT_SUCCESS) break;

      // Check whether we are still within the requested time interval.
      if (photon.time > parameters.t0 + parameters.timespan) break;

      // Check whether the photon is inside the FOV.
      // First determine telescope pointing direction at the current time.
      telescope.nz = getTelescopePointing(ac, photon.time, &status);
      if (EXIT_SUCCESS!=status) break;

      // Compare the photon direction to the unit vector specifiing the 
      // direction of the telescope axis:
      if (check_fov(&photon.direction, &telescope.nz, fov_min_align)==0) {
	// Photon is inside the FOV!
	
	// Determine telescope data like pointing direction (attitude) etc.
	// The telescope coordinate system consists of an x-, y-, and z-axis.
	// The z-axis is perpendicular to the detector plane and pointing along
	// the telescope direction. The x-axis is aligned along the detector 
	// x-direction, which is identical to the detector RAWX/COLUMN.
	// The y-axis is pointing along the y-direction of the detector,
	// which is also referred to as RAWY/ROW.

	// Determine the vector nx: perpendicular to telescope x-axis
	// and in the direction of the satellite motion.
	telescope.nx = 
	  normalize_vector(interpolate_vec(ac->entry[ac->current_entry].nx, 
					   ac->entry[ac->current_entry].time, 
					   ac->entry[ac->current_entry+1].nx, 
					   ac->entry[ac->current_entry+1].time, 
					   photon.time));
	
	// Remove the component along the vertical direction nz 
	// (nx must be perpendicular to nz!):
	scp = scalar_product(&telescope.nz, &telescope.nx);
	telescope.nx.x -= scp*telescope.nz.x;
	telescope.nx.y -= scp*telescope.nz.y;
	telescope.nx.z -= scp*telescope.nz.z;
	telescope.nx = normalize_vector(telescope.nx);


	// The third axis of the coordinate system ny is perpendicular 
	// to the telescope axis nz and nx:
	telescope.ny=normalize_vector(vector_product(telescope.nz, telescope.nx));
	
	// Determine the photon impact position on the detector (in [m]):

	// Convolution with PSF:
	// Function returns 0, if the photon does not fall on the detector. 
	// If it hits the detector, the return value is 1.
	if (get_psf_pos(&position, photon, telescope, det->focal_length, 
			det->vignetting, det->psf)) {
	  // Check whether the photon hits the detector within the FOV. 
	  // (Due to the effects of the mirrors it might have been scattered over 
	  // the edge of the FOV, although the source is inside the FOV.)
	  if (sqrt(pow(position.x,2.)+pow(position.y,2.)) < 
	      tan(det->fov_diameter)*det->focal_length) {
	    
	    // Insert the impact position with the photon data into the 
	    // impact list:
	    fits_insert_rows(impactlistfile->fptr, 
			     impactlistfile->row++, 1, &status);
	    fits_write_col(impactlistfile->fptr, TDOUBLE, impactlistfile->ctime, 
			   impactlistfile->row, 1, 1, &photon.time, &status);
	    fits_write_col(impactlistfile->fptr, TFLOAT, impactlistfile->cenergy, 
			   impactlistfile->row, 1, 1, &photon.energy, &status);
	    fits_write_col(impactlistfile->fptr, TDOUBLE, impactlistfile->cx, 
			   impactlistfile->row, 1, 1, &position.x, &status);
	    fits_write_col(impactlistfile->fptr, TDOUBLE, impactlistfile->cy, 
			   impactlistfile->row, 1, 1, &position.y, &status);
	    impactlistfile->nrows++;
	  }
	} 
	// END get_psf_pos(...)
      } 
      // End of FOV check.
    } 
    // END of scanning LOOP over the photon list.
  } while(0); // END of the error handling loop.


  // --- cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Close the FITS files.
  destroyImpactListFile(&impactlistfile, &status);
  status += closePhotonListFile(&photonlistfile);

  free_AttitudeCatalog(ac);
  destroyGenDet(&det, &status);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully!\n\n");

  return(status);
}



////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int photon_imaging_getpar(struct Parameters* parameters)
{
  int status=EXIT_SUCCESS; // Error status.

  // Get the filename of the input photon list (FITS file)
  if ((status = PILGetFname("photonlist_filename", 
			    parameters->photonlist_filename))) {
    HD_ERROR_THROW("Error reading the filename of the photon list!\n", status);
  }
  
  // Get the filename of the XML detector description.
  else if ((status = PILGetFname("xml_filename", parameters->xml_filename))) {
    HD_ERROR_THROW("Error reading the filename of the XML detector description!\n", status);
  }

  // Get the filename of the impact list output file (FITS file)
  else if ((status = PILGetFname("impactlist_filename", parameters->impactlist_filename))) {
    HD_ERROR_THROW("Error reading the filename of the impact list output file!\n", status);
  }

  // Get the start time of the simulation
  else if ((status = PILGetReal("t0", &parameters->t0))) {
    HD_ERROR_THROW("Error reading the 't0' parameter!\n", status);
  }

  // Get the timespan for the simulation
  else if ((status = PILGetReal("timespan", &parameters->timespan))) {
    HD_ERROR_THROW("Error reading the 'timespan' parameter!\n", status);
  }

  // Read the telescope number (0 for TRoPIC, 1-7 for eROSITA).
  else if ((status = PILGetInt("telescope", &parameters->telescope))) {
    HD_ERROR_THROW("Error reading the telescope number!\n", status);
  }
  if (EXIT_SUCCESS!=status) return(status);

  // Get the name of the FITS template directory.
  // First try to read it from the environment variable.
  // If the variable does not exist, read it from the PIL.
  char* buffer;
  if (NULL!=(buffer=getenv("SIXT_FITS_TEMPLATES"))) {
    strcpy(parameters->impactlist_template, buffer);
  } else {
    if ((status = PILGetFname("fits_templates", parameters->impactlist_template))) {
      HD_ERROR_THROW("Error reading the path of the FITS templates!\n", status);      
    }
  }
  // Set the impact list template file:
  strcat(parameters->impactlist_template, "/impactlist.tpl");

  return(status);
}



