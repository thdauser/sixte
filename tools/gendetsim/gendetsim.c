#include "gendetsim.h"


////////////////////////////////////
/** Main procedure. */
int gendetsim_main() {

  // Containing all programm parameters read by PIL
  struct Parameters parameters; 
  // Detector data structure (containing the pixel array, its width, ...).
  GenDet* det=NULL;
  // Input impact list.
  ImpactListFile* ilf=NULL;

  // Output event list file.
  EventListFile* elf=NULL;

  int status=EXIT_SUCCESS; // Error status.

  // Register HEATOOL:
  set_toolname("gendetsim");
  set_toolversion("0.01");

  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    headas_chat(3, "initialization ...\n");

    // Read parameters using PIL library.
    status=getpar(&parameters);
    CHECK_STATUS_BREAK(status);

    // Initialize HEADAS random number generator.
    HDmtInit(parameters.random_seed);

    // Open the FITS file with the input impact list:
    ilf=openImpactListFile(parameters.impactlist_filename,
			   READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Initialize the detector data structure.
    det=newGenDet(parameters.xml_filename, &status);
    CHECK_STATUS_BREAK(status);

    // Open the output event file.
    elf=openNewEventListFile(parameters.eventlist_filename,
			     parameters.eventlist_template,
			     &status);
    CHECK_STATUS_BREAK(status);

    // --- END of Initialization ---


    // --- Beginning of Detection Process ---

    headas_chat(3, "start detection process ...\n");

    phdetGenDet(det, ilf, elf, parameters.t0, parameters.exposure, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.

  // --- END of Detection process ---


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Destroy the detector data structure.
  destroyGenDet(&det, &status);

  // Close the event list FITS file.
  freeEventListFile(&elf, &status);

  // Close the impact list FITS file.
  freeImpactListFile(&ilf, &status);

  if (status == EXIT_SUCCESS) headas_chat(3, "finished successfully\n\n");
  return(status);
}



////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int getpar(struct Parameters* const parameters)
{
  int status=EXIT_SUCCESS; // Error status

  // Get the name of the input impact list file (FITS file).
  if ((status = PILGetFname("impactlist_filename", parameters->impactlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the input impact list file!\n", status);
    return(status);
  }

  // Get the name of the output event list file (FITS file).
  if ((status = PILGetFname("eventlist_filename", parameters->eventlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the output event list file!\n", status);
    return(status);
  }

  // Get the name of the detector XML description file (FITS file).
  if ((status = PILGetFname("xml_filename", parameters->xml_filename))) {
    HD_ERROR_THROW("Error reading the name of the detector definition XML file!\n", status);
    return(status);
  }

  // Get the start time t0.
  if ((status = PILGetReal("t0", &parameters->t0))) {
    HD_ERROR_THROW("Error reading the start time t0!\n", status);
    return(status);
  }

  // Get the time interval.
  if ((status = PILGetReal("exposure", &parameters->exposure))) {
    HD_ERROR_THROW("Error reading the simulated exposure!\n", status);
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
    strcpy(parameters->eventlist_template, buffer);
  } else {
    if ((status = PILGetFname("fits_templates", parameters->eventlist_template))) {
      HD_ERROR_THROW("Error reading the path of the FITS templates!\n", status);      
      return(status);
    }
  }
  // Set the impact list template file:
  strcat(parameters->eventlist_template, "/eventlist.tpl");

  return(status);
}


