#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include "gendetsim.h"


////////////////////////////////////
/** Main procedure. */
int gendetsim_main() {

  // Containing all programm parameters read by PIL
  struct Parameters parameters; 
  // Detector data structure (containing the pixel array, its width, ...).
  GenDet* det=NULL;

  ImpactListFile* impactlistfile=NULL;

  int status=EXIT_SUCCESS; // Error status.


  // Register HEATOOL:
  set_toolname("gendetsim");
  set_toolversion("0.01");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    headas_chat(3, "initialization ...\n");

    // Initialize HEADAS random number generator.
    HDmtInit(SIXT_HD_RANDOM_SEED);

    // Read parameters using PIL library:
    if ((status=getpar(&parameters))) break;

    // Open the FITS file with the input impact list:
    impactlistfile = openImpactListFile(parameters.impact_filename,
					READONLY, &status);
    if (EXIT_SUCCESS!=status) break;

    // Initialize the detector data structure.
    det = newGenDet(parameters.xml_filename, &status);
    if (EXIT_SUCCESS!=status) break;

    // Set the output event file.
    GenDetSetEventFile(det, parameters.event_filename, &status);
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
      if ((impact.time<parameters.t0)||(impact.time>parameters.t0+parameters.timespan)) 
	continue;

      addGenDetPhotonImpact(det, &impact, &status);
      if (EXIT_SUCCESS!=status) break;

    };
    if (EXIT_SUCCESS!=status) break;
    // END of loop over all impacts in the FITS file.
    
    // Finalize the GenDet. Perform the time-triggered operations 
    // without adding any new charges.
    operateGenDetClock(det, parameters.t0+parameters.timespan, &status);
    if (EXIT_SUCCESS!=status) break;

  } while(0); // END of the error handling loop.

  // --- END of Detection process ---


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Destroy the detector data structure.
  destroyGenDet(&det, &status);
  
  if (status == EXIT_SUCCESS) headas_chat(3, "finished successfully\n\n");
  return(status);
}



////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int getpar(struct Parameters* const parameters)
{
  int status=EXIT_SUCCESS; // Error status

  // Get the name of the input impact list file (FITS file).
  if ((status = PILGetFname("impact_filename", parameters->impact_filename))) {
    HD_ERROR_THROW("Error reading the name of the input impact list file!\n", status);
  }

  // Get the name of the output event list file (FITS file).
  else if ((status = PILGetFname("event_filename", parameters->event_filename))) {
    HD_ERROR_THROW("Error reading the name of the output event list file!\n", status);
  }

  // Get the name of the detector XML description file (FITS file).
  else if ((status = PILGetFname("xml_filename", parameters->xml_filename))) {
    HD_ERROR_THROW("Error reading the name of the detector definition XML file!\n", status);
  }

  // Get the start time t0.
  else if ((status = PILGetReal("t0", &parameters->t0))) {
    HD_ERROR_THROW("Error reading the start time t0!\n", status);
  }

  // Get the time interval.
  else if ((status = PILGetReal("timespan", &parameters->timespan))) {
    HD_ERROR_THROW("Error reading the simulated timespan!\n", status);
  }

  return(status);
}


