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

  ImpactListFile* impactlistfile=NULL;

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
    HDmtInit(SIXT_HD_RANDOM_SEED);


    // Open the impact list FITS file.
    impactlistfile = openImpactListFile(parameters.impactlist_filename, 
					READONLY, &status);
    if (EXIT_SUCCESS!=status) break;


    // Detector settings.
    // Store the settings for the XMSDetector in the corresponding data structure
    // and call the initialization routine in order to allocate memory etc.
    struct XMSDetectorParameters hdparameters = {
      .pixels_inner = { .xwidth=parameters.width_inner,
			.ywidth=parameters.width_inner,
			.xpixelwidth = parameters.pixelwidth_inner,
			.ypixelwidth = parameters.pixelwidth_inner },
      .pixels_outer = { .xwidth=parameters.width_outer,
			.ywidth=parameters.width_outer,
			.xpixelwidth = parameters.pixelwidth_outer,
			.ypixelwidth = parameters.pixelwidth_outer },
      .generic_inner = { .ccsigma = parameters.ccsigma_inner, 
			 .pha_threshold = parameters.pha_threshold_inner,
			 .energy_threshold = parameters.energy_threshold_inner,
			 .rmf_filename = parameters.rmf_filename_inner /* String address!! */ },
      .generic_outer = { .ccsigma = parameters.ccsigma_outer, 
			 .pha_threshold = parameters.pha_threshold_outer,
			 .energy_threshold = parameters.energy_threshold_outer,
			 .rmf_filename = parameters.rmf_filename_outer /* String address!! */ },
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
    while ((EXIT_SUCCESS==status)&&(0==ImpactListFile_EOF(impactlistfile))) {

      getNextImpactFromFile(impactlistfile, &impact, &status);
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

  // Close the FITS files.
  destroyImpactListFile(&impactlistfile, &status);

  // Release memory of detector.
  status+=cleanupXMSDetector(&detector);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully\n\n");
  return(status);
}




////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
static int getpar(struct Parameters* parameters)
{
  int pha_threshold;
  int status=EXIT_SUCCESS; // error status

  // Get the name of the impact list file (FITS file):
  if ((status = PILGetFname("impactlist_filename", parameters->impactlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the impact list file!\n", status);
  }


  // INNER TES pixel array:
  // Get the detector width (number of pixels):
  else if ((status = PILGetInt("det_width_inner", &parameters->width_inner))) {
    HD_ERROR_THROW("Error reading the width of the inner TES pixel array!\n", status);
  }

  // [m]
  else if ((status = PILGetReal("pixelwidth_inner", &parameters->pixelwidth_inner))) {
    HD_ERROR_THROW("Error reading the width of the inner detector pixels!\n", status);
  }

  // [m]
  else if ((status = PILGetReal("ccsigma_inner", &parameters->ccsigma_inner))) {
    HD_ERROR_THROW("Error reading the charge cloud sigma for the inner pixels!\n", status);
  }

  
  // Read the detector thresholds (either integer PHA or float energy):
  else if ((status = PILGetInt("pha_threshold_inner", &pha_threshold))) {
    HD_ERROR_THROW("Error: could not determine detector PHA threshold for inner pixels!\n", 
		   status);
    return(status);
  } else {
    parameters->pha_threshold_inner = (long)pha_threshold;
  }
  if (parameters->pha_threshold_inner==-1) {
    if ((status = PILGetReal4("energy_threshold_inner", &parameters->energy_threshold_inner))) {
      HD_ERROR_THROW("Error: could not determine detector energy threshold for inner pixesl!\n", 
		     status);
      return(status);
    }
  } else {
    parameters->energy_threshold_inner=0.;
  }

  // Get the name of the detector redistribution file (FITS file)
  if ((status = PILGetFname("rmf_filename_inner", parameters->rmf_filename_inner))) {
    HD_ERROR_THROW("Error reading the name of the detector redistribution matrix"
		   "file (RMF) for the inner TES pixel array!\n", status);
  }


  // OUTER TES pixel array:
  // Get the detector width (number of pixels):
  else if ((status = PILGetInt("det_width_outer", &parameters->width_outer))) {
    HD_ERROR_THROW("Error reading the width of the outer TES pixel array!\n", status);
  }

  // [m]
  else if ((status = PILGetReal("pixelwidth_outer", &parameters->pixelwidth_outer))) {
    HD_ERROR_THROW("Error reading the width of the outer detector pixels!\n", status);
  }

  // [m]
  else if ((status = PILGetReal("ccsigma_outer", &parameters->ccsigma_outer))) {
    HD_ERROR_THROW("Error reading the charge cloud sigma for the outer pixels!\n", status);
  }

  
  // Read the detector thresholds (either integer PHA or float energy):
  else if ((status = PILGetInt("pha_threshold_outer", &pha_threshold))) {
    HD_ERROR_THROW("Error: could not determine detector PHA threshold for outer pixels!\n", 
		   status);
    return(status);
  } else {
    parameters->pha_threshold_outer = (long)pha_threshold;
  }
  if (parameters->pha_threshold_outer==-1) {
    if ((status = PILGetReal4("energy_threshold_outer", &parameters->energy_threshold_outer))) {
      HD_ERROR_THROW("Error: could not determine detector energy threshold for outer pixesl!\n", 
		     status);
      return(status);
    }
  } else {
    parameters->energy_threshold_outer=0.;
  }

  // Get the name of the detector redistribution file (FITS file)
  if ((status = PILGetFname("rmf_filename_outer", parameters->rmf_filename_outer))) {
    HD_ERROR_THROW("Error reading the name of the detector redistribution matrix"
		   "file (RMF) for the outer TES pixel array!\n", status);
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


