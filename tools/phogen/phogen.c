#include "phogen.h"


int phogen_main() 
{
  // Program parameters.
  struct Parameters parameters;
  
  // Error status.
  int status = EXIT_SUCCESS;  

  // Register HEATOOL.
  set_toolname("phogen");
  set_toolversion("0.02");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    status=phogen_getpar(&parameters);
    CHECK_STATUS_BREAK(status);

    // Initialize HEADAS random number generator.
    HDmtInit(SIXT_HD_RANDOM_SEED);
    
    // --- End of Initialization ---


    // --- Photon Generation Process ---

    // Start the actual photon generation (after loading required data):
    headas_chat(5, "start photon generation process ...\n");

    photon_generation(parameters.xml_filename,
		      parameters.attitude_filename,
		      parameters.simput_filename,
		      parameters.photonlist_filename,
		      parameters.t0, 
		      parameters.t0+parameters.timespan,
		      &status);
    CHECK_STATUS_BREAK(status);

    // --- End of photon generation ---

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(3, "\ncleaning up ...\n");

  // Release HEADAS random number generator:
  HDmtFree();

  if (status==EXIT_SUCCESS) headas_chat(0, "finished successfully!\n\n");
  return(status);
}



int phogen_getpar(struct Parameters* parameters)
{
  char msg[MAXMSG];           // Error message buffer.
  int status = EXIT_SUCCESS;  // Error status flag.

  // Get the filename of the detector XML definition file.
  if ((status = PILGetFname("xml_filename", parameters->xml_filename))) {
    HD_ERROR_THROW("Error reading the filename of the detector " 
		   "XML definition file!\n", status);
  }

  // Get the filename of the Attitude file (FITS file).
  else if ((status = PILGetFname("attitude_filename", parameters->attitude_filename))) {
    HD_ERROR_THROW("Error reading the filename of the attitude file!\n", status);
  }

  // Get the start time of the photon generation
  else if ((status = PILGetReal("t0", &parameters->t0))) {
    HD_ERROR_THROW("Error reading the 't0' parameter!\n", status);
  }

  // Get the timespan for the photon generation
  else if ((status = PILGetReal("timespan", &parameters->timespan))) {
    HD_ERROR_THROW("Error reading the 'timespan' parameter!\n", status);
  }

  // Determine the name of the file that contains the input sources (either
  // a point source catalog, source images, or a FITS grouping extension listing
  // several input files).
  else if ((status = PILGetFname("simput_filename", parameters->simput_filename))) {
    HD_ERROR_THROW("Error reading the filename of the input sources!\n", status);
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
  if (EXIT_SUCCESS!=status) return(status);
  // Set the photon list template file:
  strcat(parameters->photonlist_template, "/photonlist.tpl");

  return(status);
}


