#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include "comadet.h"


////////////////////////////////////
/** Main procedure. */
int comadet_main() {
  struct Parameters parameters;
  
  ImpactListFile* impactlistfile=NULL;
  CoMaDetector* detector=NULL;

  int status=EXIT_SUCCESS; // Error status.
  //  char msg[MAXMSG];


  // Register HEATOOL:
  set_toolname("comadet");
  set_toolversion("0.01");


  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // Read the program parameters using the PIL library.
    if ((status=comadet_getpar(&parameters))) break;
    
    // Open the impact list FITS file.
    impactlistfile = openImpactListFile(parameters.impactlist_filename,
					READONLY, &status);
    if (EXIT_SUCCESS!=status) break;
    
    // DETECTOR setup.
    struct CoMaDetectorParameters cdp = {
      .pixels = 
      { .xwidth = parameters.width,
	.ywidth = parameters.width,
	.xpixelwidth = parameters.pixelwidth,
	.ypixelwidth = parameters.pixelwidth 
      },
      .eventfile_filename = parameters.eventlist_filename /* String address!! */,
      .eventfile_template = parameters.eventlist_template
    };
    detector=getCoMaDetector(&cdp, &status);
    if(EXIT_SUCCESS!=status) break;
    // END of DETECTOR CONFIGURATION SETUP    

    // --- END of Initialization ---


    // --- Beginning of Imaging Process ---

    // Beginning of actual detector simulation (after loading required data):
    headas_chat(5, "start detection process ...\n");
    Impact impact;

    // Loop over all impacts in the FITS file.
    while ((EXIT_SUCCESS==status)&&(0==ImpactListFile_EOF(impactlistfile))) {

      status=getNextImpactListFileRow(impactlistfile, &impact);
      if(EXIT_SUCCESS!=status) break;

      // Call the photon detection routine of the Coded Mask Detector.
      status=addImpact2CoMaDetector(detector, &impact);
      if (EXIT_SUCCESS!=status) break;

    } // END of scanning the impact list.
    if (EXIT_SUCCESS!=status) break;

  } while(0);  // END of the error handling loop.


  // --- cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Free the Coded Mask Detector.
  freeCoMaDetector(detector);

  // Close the FITS files.
  destroyImpactListFile(&impactlistfile, &status);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully!\n\n");
  return(status);
}



int comadet_getpar(struct Parameters* parameters)
{
  int status=EXIT_SUCCESS; // Error status.

  // Get the filename of the impact list file (FITS input file).
  if ((status = PILGetFname("impactlist_filename", 
			    parameters->impactlist_filename))) {
    HD_ERROR_THROW("Error reading the filename of the impact list intput "
		   "file!\n", status);
  }

  // Get the filename of the event list file (FITS output file).
  else if ((status = PILGetFname("eventlist_filename", 
				 parameters->eventlist_filename))) {
    HD_ERROR_THROW("Error reading the filename of the event list output "
		   "file!\n", status);
  }

  // Read the width of the detector in [pixel].
  else if ((status = PILGetInt("width", &parameters->width))) {
    HD_ERROR_THROW("Error reading the detector width!\n", status);
  }

  // Read the width of one detector pixel in [m].
  else if ((status = PILGetReal("pixelwidth", &parameters->pixelwidth))) {
    HD_ERROR_THROW("Error reading the width detector pixels!\n", status);
  }
  if (EXIT_SUCCESS!=status) return(status);

  // Get the name of the FITS template directory.
  // First try to read it from the environment variable.
  // If the variable does not exist, read it from the PIL.
  char* buffer;
  if (NULL!=(buffer=getenv("SIXT_FITS_TEMPLATES"))) {
    strcpy(parameters->eventlist_template, buffer);
  } else {
    if ((status = PILGetFname("fits_templates", 
			      parameters->eventlist_template))) {
      HD_ERROR_THROW("Error reading the path of the FITS templates!\n", status);
      
    }
  }
  // Set the impact list template file:
  strcat(parameters->eventlist_template, "/coma.eventlist.tpl");

  return(status);
}



