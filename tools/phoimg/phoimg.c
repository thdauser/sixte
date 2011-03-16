#include "phoimg.h"


////////////////////////////////////
/** Main procedure. */
int phoimg_main() {
  struct Parameters parameters;

  AttitudeCatalog* ac=NULL;

  PhotonListFile* photonlistfile=NULL;
  ImpactListFile* impactlistfile=NULL;

  // Detector data structure including telescope information like the PSF,
  // vignetting function, focal length, and FOV diameter.
  GenDet* det=NULL;

  int status=EXIT_SUCCESS; // Error status


  // Register HEATOOL:
  set_toolname("phoimg");
  set_toolversion("0.01");


  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // Read parameters using PIL library.
    if ((status=phoimg_getpar(&parameters))) break;

    // Initialize the detector data structure.
    det = newGenDet(parameters.xml_filename, &status);
    CHECK_STATUS_BREAK(status);
    
    // Initialize HEADAS random number generator.
    HDmtInit(parameters.random_seed);

    // Open the FITS file with the input photon list:
    photonlistfile=openPhotonListFile(parameters.photonlist_filename, 
				      READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Open the attitude file specified in the header keywords of the photon list.
    char comment[MAXMSG]; // String buffer.
    if (fits_read_key(photonlistfile->fptr, TSTRING, "ATTITUDE", 
		      &parameters.attitude_filename, 
		      comment, &status)) break;
    if (0<strlen(parameters.attitude_filename)) {
      if (NULL==(ac=loadAttitudeCatalog(parameters.attitude_filename,
					parameters.t0, parameters.exposure, 
					&status))) break;
    } else {
      status=EXIT_FAILURE;
      HD_ERROR_THROW("No information about attitude file!\n", status);
      break;
    }

    // Create a new FITS file for the output of the impact list.
    impactlistfile = openNewImpactListFile(parameters.impactlist_filename, 
					   parameters.impactlist_template,
					   &status);
    CHECK_STATUS_BREAK(status);

    // Add attitude filename.
    if (fits_update_key(impactlistfile->fptr, TSTRING, "ATTITUDE", parameters.attitude_filename,
		       "name of the attitude FITS file", &status)) break;
    
    // --- END of Initialization ---


    // --- Beginning of Imaging Process ---

    // Beginning of actual simulation (after loading required data):
    headas_chat(3, "start imaging process ...\n");

    phimg(det, ac, photonlistfile, impactlistfile,
	  parameters.t0, parameters.exposure, &status);
    CHECK_STATUS_BREAK(status);

    // --- END of imaging process ---

  } while(0); // END of the error handling loop.


  // --- cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Close the FITS files.
  freeImpactListFile(&impactlistfile, &status);
  freePhotonListFile(&photonlistfile, &status);

  freeAttitudeCatalog(&ac);
  destroyGenDet(&det, &status);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully!\n\n");

  return(status);
}



////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int phoimg_getpar(struct Parameters* parameters)
{
  int status=EXIT_SUCCESS; // Error status.

  // Get the filename of the input photon list (FITS file)
  if ((status = PILGetFname("photonlist_filename", 
			    parameters->photonlist_filename))) {
    HD_ERROR_THROW("Error reading the filename of the photon list!\n", status);
    return(status);
  }
  
  // Get the filename of the XML detector description.
  if ((status = PILGetFname("xml_filename", parameters->xml_filename))) {
    HD_ERROR_THROW("Error reading the filename of the XML detector description!\n", status);
    return(status);
  }

  // Get the filename of the impact list output file (FITS file)
  if ((status = PILGetFname("impactlist_filename", parameters->impactlist_filename))) {
    HD_ERROR_THROW("Error reading the filename of the impact list output file!\n", status);
    return(status);
  }

  // Get the start time of the simulation
  if ((status = PILGetReal("t0", &parameters->t0))) {
    HD_ERROR_THROW("Error reading the 't0' parameter!\n", status);
    return(status);
  }

  // Get the timespan for the simulation
  if ((status = PILGetReal("exposure", &parameters->exposure))) {
    HD_ERROR_THROW("Error reading the 'exposure' parameter!\n", status);
    return(status);
  }

  // Get the seed for the random number generator.
  if ((status = PILGetInt("random_seed", &parameters->random_seed))) {
    HD_ERROR_THROW("Error reading the seed for the random "
		   "number generator!\n", status);
    return(status);
  }

  // Get the name of the FITS template directory.
  // First try to read it from the environment variable.
  // If the variable does not exist, read it from the PIL.
  char* buffer;
  if (NULL!=(buffer=getenv("SIXT_FITS_TEMPLATES"))) {
    strcpy(parameters->impactlist_template, buffer);
  } else {
    if ((status = PILGetFname("fits_templates", parameters->impactlist_template))) {
      HD_ERROR_THROW("Error reading the path of the FITS templates!\n", status);      
      return(status);
    }
  }
  // Set the impact list template file:
  strcat(parameters->impactlist_template, "/impactlist.tpl");

  return(status);
}



