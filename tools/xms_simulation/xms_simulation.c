#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include "xms_simulation.h"


////////////////////////////////////
/** Main procedure. */
int xms_simulation_main() {

  // Containing all programm parameters read by PIL
  struct Parameters parameters; 

  // XMSDetector data structure.
  // This is the main object for the simulation of the XMS. It contains all
  // important parameters of the detector and is used to create an event file
  // from the incident photon events.
  // Before the first usage it has to be initialized.
  XMSDetector detector;

  struct ImpactlistFile impactlistfile;

  int status=EXIT_SUCCESS; // Error status.


  // Register HEATOOL:
  set_toolname("xms_simulation");
  set_toolversion("0.01");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    // Read parameters using PIL library:
    if ((status=getpar(&parameters))) break;
    
    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);


    // Open the impact list FITS file.
    if(EXIT_SUCCESS!=(status=impactlist_openFile(&impactlistfile, 
						 parameters.impactlist_filename, 
						 READONLY))) break;

    // Detector settings.
    // Store the settings for the XMSDetector in the corresponding data structure
    // and call the initialization routine in order to allocate memory etc.
    struct XMSDetectorParameters hdparameters = {
      .pixels = { .xwidth=parameters.width,
		  .ywidth=parameters.width,
		  .xpixelwidth = parameters.pixelwidth,
		  .ypixelwidth = parameters.pixelwidth },
      .generic = { .ccsigma = parameters.ccsigma, 
		   .pha_threshold = parameters.pha_threshold,
		   .energy_threshold = parameters.energy_threshold,
		   .rmf_filename = parameters.rmf_filename /* String address!! */ },
      .eventlist_filename = parameters.eventlist_filename /* String address!! */,
      .eventlist_template = parameters.eventlist_template /* String address!! */
    };
    if(EXIT_SUCCESS!=(status=initXMSDetector(&detector, &hdparameters))) break;
    // END of DETECTOR CONFIGURATION SETUP    

    // --- END of Initialization ---


    // --- Beginning of Detection Process ---

    headas_chat(5, "start detection process ...\n");
    Impact impact; // Buffer to store the impacts read from the FITS file.

    // Loop over all impacts in the FITS file.
    while ((EXIT_SUCCESS==status)&&(0==impactlist_EOF(&impactlistfile))) {

      status=impactlist_getNextRow(&impactlistfile, &impact);
      if (EXIT_SUCCESS!=status) break;

      // Check whether the event lies in the specified time interval:
      if ((impact.time<parameters.t0)||(impact.time>parameters.t0+parameters.timespan)) continue;

      // Call the photon detection routine that generates the right charge
      // and stores it in the detector pixels.
      // Before generating and adding the charge to the detector the routine also 
      // checks, whether the integration time is exceeded and performs the readout 
      // in that case. 
      status=addImpact2XMSDetector(&detector, &impact);
      if (EXIT_SUCCESS!=status) break;

    } // END of scanning the impact list.
    
  } while(0); // END of the error handling loop.

  // --- END of Detection process ---


  // --- Cleaning up ---
  headas_chat(5, "\ncleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  if (NULL!=impactlistfile.fptr) fits_close_file(impactlistfile.fptr, &status);

  // Release memory of detector.
  status+=cleanupXMSDetector(&detector);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully\n\n");
  return(status);
}




////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
static int getpar(struct Parameters* parameters)
{
  int status=EXIT_SUCCESS; // error status

  // Get the name of the impact list file (FITS file):
  if ((status = PILGetFname("impactlist_filename", parameters->impactlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the impact list file!\n", status);
  }

  // Get the detector width (number of pixels):
  else if ((status = PILGetInt("det_width", &parameters->width))) {
    HD_ERROR_THROW("Error reading the width of the detector!\n", status);
  }

  // [m]
  else if ((status = PILGetReal("pixelwidth", &parameters->pixelwidth))) {
    HD_ERROR_THROW("Error reading the width of the detector pixels!\n", status);
  }

  // [m]
  else if ((status = PILGetReal("ccsigma", &parameters->ccsigma))) {
    HD_ERROR_THROW("Error reading the charge cloud sigma!\n", status);
  }
  if (status) return(status);

  
  // Read the detector thresholds (either integer PHA or float energy):
  int pha_threshold;
  if ((status = PILGetInt("pha_threshold", &pha_threshold))) {
    HD_ERROR_THROW("Error: could not determine detector PHA threshold!\n", status);
    return(status);
  } else {
    parameters->pha_threshold = (long)pha_threshold;
  }
  if (parameters->pha_threshold==-1) {
    if ((status = PILGetReal4("energy_threshold", &parameters->energy_threshold))) {
      HD_ERROR_THROW("Error: could not determine detector energy threshold!\n", status);
      return(status);
    }
  } else {
    parameters->energy_threshold=0.;
  }

  // Get the name of the detector redistribution file (FITS file)
  if ((status = PILGetFname("rmf_filename", parameters->rmf_filename))) {
    HD_ERROR_THROW("Error reading the name of the detector" 
		   "redistribution matrix file (RMF)!\n", status);
  }

  // Get the background count rate
  else if ((status = PILGetReal4("background_rate", &parameters->background_rate))) {
    HD_ERROR_THROW("Error: could not determine the detector background rate!\n", status);
  }

  // Get the name of the output event list (FITS file)
  else if ((status = PILGetFname("eventlist_filename", parameters->eventlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the event list file!\n", status);
  }

  // Get the start time of the simulation
  else if ((status = PILGetReal("t0", &parameters->t0))) {
    HD_ERROR_THROW("Error reading the 't0' parameter!\n", status);
  }

  // Get the timespan for the simulation
  else if ((status = PILGetReal("timespan", &parameters->timespan))) {
    HD_ERROR_THROW("Error reading the 'timespan' parameter!\n", status);
  }
  
  // Get the name of the FITS template directory.
  // First try to read it from the environment variable.
  // If the variable does not exist, read it from the PIL.
  else { 
    char* buffer;
    if (NULL!=(buffer=getenv("SIXT_FITS_TEMPLATES"))) {
      strcpy(parameters->eventlist_template, buffer);
    } else {
      if ((status = PILGetFname("fits_templates", parameters->eventlist_template))) {
	HD_ERROR_THROW("Error reading the path of the FITS templates!\n", status);
      }
    }
  }

  // Set the event list template file:
  strcat(parameters->eventlist_template, "/xms.eventlist.tpl");

  return(status);
}


