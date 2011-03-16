#include "gendetsim.h"


////////////////////////////////////
/** Main procedure. */
int gendetsim_main() {

  // Containing all programm parameters read by PIL
  struct Parameters parameters; 
  // Detector data structure (containing the pixel array, its width, ...).
  GenDet* det=NULL;
  // Input impact list.
  ImpactListFile* impactlistfile=NULL;
  // Total number of detected photons. Only the number of
  // photons absorbed by valid pixels inside the detector is
  // counted. Split events created by one photon are counted only
  // once.
  long n_detected_photons=0;

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
    if ((status=getpar(&parameters))) break;

    // Initialize HEADAS random number generator.
    HDmtInit(parameters.random_seed);

    // Open the FITS file with the input impact list:
    impactlistfile = openImpactListFile(parameters.impactlist_filename,
					READONLY, &status);
    if (EXIT_SUCCESS!=status) break;

    // Initialize the detector data structure.
    det = newGenDet(parameters.xml_filename, &status);
    if (EXIT_SUCCESS!=status) break;

    // Open the output event file.
    elf = openNewEventListFile(parameters.eventlist_filename,
			       parameters.eventlist_template,
			       &status);
    if (EXIT_SUCCESS!=status) break;

    // --- END of Initialization ---


    // --- Beginning of Detection Process ---

    headas_chat(3, "start detection process ...\n");

    // Loop over all impacts in the FITS file.
    Impact impact;
    while ((EXIT_SUCCESS==status)&&(0==ImpactListFile_EOF(impactlistfile))) {

      getNextImpactFromFile(impactlistfile, &impact, &status);
      if(EXIT_SUCCESS!=status) break;

      // Check whether the event lies in the specified time interval:
      if ((impact.time<parameters.t0)||(impact.time>parameters.t0+parameters.exposure)) 
	continue;

      // Add the impact to the detector array. If it is absorbed
      // by at least one valid pixel, increase the counter for
      // the number of detected photons.
      if (addGenDetPhotonImpact(det, &impact, elf, &status) > 0) {
	n_detected_photons++;
      }
      if (EXIT_SUCCESS!=status) break;

    };
    if (EXIT_SUCCESS!=status) break;
    // END of loop over all impacts in the FITS file.
    
    // Finalize the GenDet. Perform the time-triggered operations 
    // without adding any new charges.
    operateGenDetClock(det, elf, parameters.t0+parameters.exposure, &status);
    if (EXIT_SUCCESS!=status) break;

    // Store the number of simulated input photons in the FITS header
    // of the output event file.
    if (fits_update_key(elf->fptr, TLONG, "NPHOTONS", 
			&impactlistfile->nrows, "number of input photons", 
			&status)) break;
    // Store the number of detected photons in the FITS header of
    // the output event file.
    if (fits_update_key(elf->fptr, TLONG, "NDETECTD", 
			&n_detected_photons, "number of detected photons", 
			&status)) break;

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
  freeImpactListFile(&impactlistfile, &status);

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


