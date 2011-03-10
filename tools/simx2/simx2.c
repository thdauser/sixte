#include "simx2.h"


int simx2_main() 
{
  // Program parameters.
  struct Parameters par;
  
  int status = EXIT_SUCCESS; // Error status.

  // Register HEATOOL
  set_toolname("simx2");
  set_toolversion("0.01");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----

    status=simx2_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // --- End of Initialization ---


    // --- Simulation Process ---


    // --- End of simulation process ---

  } while(0); // END of ERROR HANDLING Loop.


  // --- Clean up ---
  
  headas_chat(3, "\ncleaning up ...\n");

  if (status==EXIT_SUCCESS) headas_chat(0, "finished successfully!\n\n");
  return(status);
}



int simx2_getpar(struct Parameters* const par)
{
  int status = EXIT_SUCCESS; // Error status.

  // Get the filename of the detector XML definition file.
  if ((status = PILGetFname("xml_filename", par->xml_filename))) {
    HD_ERROR_THROW("Error reading the name of the detector " 
		   "XML definition file!\n", status);
    return(status);
  }

  // Get the filename of the attitude file (FITS file).
  if ((status = PILGetString("attitude_filename", par->attitude_filename))) {
    HD_ERROR_THROW("Error reading the name of the attitude file!\n", status);
    return(status);

  } else if (0==strlen(par->attitude_filename)) {
    // If the attitude filename is empty, read pointing parameters 
    // from the command line / PIL.
    if ((status = PILGetReal4("pointing_ra", &par->pointing_ra))) {
      HD_ERROR_THROW("Error reading the right ascension of the telescope "
		     "pointing direction!\n", status);
      return(status);
    }
    if ((status = PILGetReal4("pointing_dec", &par->pointing_dec))) {
      HD_ERROR_THROW("Error reading the declination of the telescope "
		     "pointing direction!\n", status);
      return(status);
    }
  }

  // Determine the name of the file that contains the input source catalog.
  // The file must have the SIMPUT format.
  if ((status = PILGetString("simput_filename", par->simput_filename))) {
    HD_ERROR_THROW("Error reading the filename of the input source "
		   "catalog (SIMPUT)!\n", status);
  }
  // TODO If the source catalog filename is empty, read source parameters
  // from the command line / PIL.

  // Get the start time of the simulation.
  if ((status = PILGetReal("t0", &par->t0))) {
    HD_ERROR_THROW("Error reading the start time t_0!\n", status);
    return(status);
  }

  // Get the exposure time for the simulation.
  if ((status = PILGetReal("exposure", &par->exposure))) {
    HD_ERROR_THROW("Error reading the exposure time!\n", status);
    return(status);
  }

  // Get the filename of the output event list file.
  if ((status = PILGetFname("eventlist_filename", par->eventlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the output " 
		   "event list file!\n", status);
    return(status);
  }

  // Get the filename of the output event list file.
  if ((status = PILGetInt("random_seed", &par->random_seed))) {
    HD_ERROR_THROW("Error reading the seed for the random "
		   "number generator!\n", status);
    return(status);
  }

  // Get the name of the FITS template directory from the 
  // environment variable.
  char* buffer;
  if (NULL!=(buffer=getenv("SIXT_FITS_TEMPLATES"))) {
    strcpy(par->fits_templates, buffer);
  } else {
    // Could not read the environment variable.
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error reading the environment variable containing "
		   "the location of the FITS templates!\n", status);
    return(status);
  }

  return(status);
}


