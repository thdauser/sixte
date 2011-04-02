#include "phogen.h"


int phogen_main() 
{
  // Program parameters.
  struct Parameters parameters;

  // Detector setup.
  GenDet* det=NULL;

  // Attitude.
  AttitudeCatalog* ac=NULL;

  // Catalog of input X-ray sources.
  SourceCatalog* srccat=NULL;

  // Photon list file.
  PhotonListFile* plf=NULL;
  
  // Error status.
  int status = EXIT_SUCCESS;  

  // Register HEATOOL.
  set_toolname("phogen");
  set_toolversion("0.02");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    // Read the parameters using PIL.
    status=phogen_getpar(&parameters);
    CHECK_STATUS_BREAK(status);

    // Initialize HEADAS random number generator.
    HDmtInit(parameters.random_seed);
    
    // Load the detector configuration.
    det=newGenDet(parameters.xml_filename, &status);
    CHECK_STATUS_BREAK(status);

    // Load the attitude from the given file.
    ac=loadAttitudeCatalog(parameters.attitude_filename,
			   parameters.t0, parameters.exposure, &status);
    CHECK_STATUS_BREAK(status);

    // Load the SIMPUT X-ray source catalog.
    srccat=loadSourceCatalog(parameters.simput_filename, det, &status);
    CHECK_STATUS_BREAK(status);

    // Remove the old photon list file.    
    remove(parameters.photonlist_filename);
    
    // Open the output photon list file.
    plf=openNewPhotonListFile(parameters.photonlist_filename, 
			      parameters.photonlist_template, &status);
    CHECK_STATUS_BREAK(status);
    // Set the attitude filename in the photon list (obsolete).
    char buffer[MAXMSG];
    strcpy(buffer, parameters.attitude_filename);
    fits_update_key(plf->fptr, TSTRING, "ATTITUDE", buffer,
		    "attitude file", &status);
    CHECK_STATUS_BREAK(status);

    // --- End of Initialization ---


    // --- Photon Generation Process ---

    // Start the actual photon generation (after loading required data):
    headas_chat(3, "start photon generation process ...\n");

    phgen(det, ac, srccat, plf, parameters.t0, parameters.exposure, 
	  &status);
    CHECK_STATUS_BREAK(status);
 
    // --- End of photon generation ---

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(3, "\ncleaning up ...\n");

  // Release memory.
  freePhotonListFile(&plf, &status);
  freeSourceCatalog(&srccat);
  freeAttitudeCatalog(&ac);
  destroyGenDet(&det, &status);

  // Release HEADAS random number generator:
  HDmtFree();

  if (status==EXIT_SUCCESS) headas_chat(0, "finished successfully!\n\n");
  return(status);
}



int phogen_getpar(struct Parameters* parameters)
{
  int status = EXIT_SUCCESS; // Error status.

  // Get the filename of the detector XML definition file.
  if ((status = PILGetFname("xml_filename", parameters->xml_filename))) {
    HD_ERROR_THROW("Error reading the filename of the detector " 
		   "XML definition file!\n", status);
    return(status);
  }

  // Get the filename of the Attitude file (FITS file).
  if ((status = PILGetFname("attitude_filename", parameters->attitude_filename))) {
    HD_ERROR_THROW("Error reading the filename of the attitude file!\n", status);
    return(status);
  }

  // Get the start time of the photon generation
  if ((status = PILGetReal("t0", &parameters->t0))) {
    HD_ERROR_THROW("Error reading the 't0' parameter!\n", status);
    return(status);
  }

  // Get the timespan for the photon generation
  if ((status = PILGetReal("exposure", &parameters->exposure))) {
    HD_ERROR_THROW("Error reading the 'exposure' parameter!\n", status);
    return(status);
  }

  // Determine the name of the file that contains the input sources (either
  // a point source catalog, source images, or a FITS grouping extension listing
  // several input files).
  if ((status = PILGetFname("simput_filename", parameters->simput_filename))) {
    HD_ERROR_THROW("Error reading the filename of the input sources!\n", status);
    return(status);
  }

  // Get the filename of the Photon-List file (FITS file):
  if ((status = PILGetFname("photonlist_filename", 
			    parameters->photonlist_filename))) {
    HD_ERROR_THROW("Error reading the filename of the output file for "
		   "the photon list!\n", status);
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
    strcpy(parameters->photonlist_template, buffer);
  } else {
    if ((status = PILGetFname("fits_templates", parameters->photonlist_template))) {
      HD_ERROR_THROW("Error reading the path of the FITS templates!\n", status);
      return(status);      
    }
  }
  // Set the photon list template file:
  strcat(parameters->photonlist_template, "/photonlist.tpl");

  return(status);
}


