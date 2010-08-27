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

  int status=EXIT_SUCCESS; // Error status.


  // Register HEATOOL:
  set_toolname("gendetsim");
  set_toolversion("0.01");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    headas_chat(5, "initialization ...\n");

    // Initialize HEADAS random number generator.
    HDmtInit(SIXT_HD_RANDOM_SEED);

    // Read parameters using PIL library:
    if ((status=getpar(&parameters))) break;

    // Initialize the detector data structure.
    det = newGenDet(parameters.xml_filename, &status);
    if (EXIT_SUCCESS!=status) break;

    // --- END of Initialization ---


    // --- Beginning of Detection Process ---

    headas_chat(5, "start detection process ...\n");

    Impact impact = { .position = { .x = 0.0144,
				    .y = 0. },
		      .time = 0.,
		      .energy = 0. };

    addGenDetPhotonImpact(det, &impact, &status);

  } while(0); // END of the error handling loop.

  // --- END of Detection process ---


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Destroy the detector data structure.
  destroyGenDet(&det, &status);
  
  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully\n\n");
  return(status);
}



////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int getpar(struct Parameters* const parameters)
{
  int status=EXIT_SUCCESS; // Error status

  // Get the name of the impact list file (FITS file)
  if ((status = PILGetFname("xml_filename", parameters->xml_filename))) {
    HD_ERROR_THROW("Error reading the name of the detector definition XML file!\n", status);
  }

  return(status);
}


