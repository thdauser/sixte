#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include "framestore_simulation.h"


////////////////////////////////////
/** Main procedure. */
int framestore_simulation_main() {
  struct Parameters parameters; // Containing all programm parameters read by PIL

  // Detector data structure (containing the pixel array, its width, ...)
  FramestoreDetector detector;

  struct ImpactlistFile impactlistfile;

  int status=EXIT_SUCCESS;   // error status


  // Register HEATOOL:
  set_toolname("framestore_simulation");
  set_toolversion("0.01");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    // Read parameters using PIL library:
    if ((status=getpar(&parameters))) break;


    // Open the impact list FITS file.
    if(EXIT_SUCCESS!=(status=impactlist_openFile(&impactlistfile, parameters.impactlist_filename, 
						 READONLY))) break;
    
    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);


    // General detector settings.
    struct FramestoreDetectorParameters fdparameters = {
      .pixels = { .xwidth = parameters.width,
		  .ywidth = parameters.width,
		  .xpixelwidth = parameters.pixelwidth,
		  .ypixelwidth = parameters.pixelwidth },
      .generic = { .ccsigma = parameters.ccsigma, 
		   .pha_threshold = parameters.pha_threshold,
		   .energy_threshold = parameters.energy_threshold,
		   .rmf_filename = parameters.rmf_filename /* String address!! */ },
      .integration_time   = parameters.integration_time,
      .t0                 = parameters.t0,
      .eventlist_filename = parameters.eventlist_filename /* String address!! */,
      .eventlist_template = parameters.eventlist_template
    };
    if(EXIT_SUCCESS!=(status=initFramestoreDetector(&detector, &fdparameters))) break;
    
    // END of DETECTOR CONFIGURATION SETUP    

    // Add important additional HEADER keywords to the event list.
    if (fits_write_key(detector.eventlist.generic.fptr, TSTRING, "ATTITUDE", 
		       impactlistfile.attitude_filename,
		       "name of the attitude FITS file", &status)) break;

    // --- END of Initialization ---


    // --- Beginning of Detection Process ---

    headas_chat(5, "start detection process ...\n");
    Impact impact;

    // Loop over all impacts in the FITS file.
    while ((EXIT_SUCCESS==status)&&(0==impactlist_EOF(&impactlistfile))) {

      status=impactlist_getNextRow(&impactlistfile, &impact);
      if(EXIT_SUCCESS!=status) break;

      // Check whether the event lies in the specified time interval:
      if ((impact.time<parameters.t0)||(impact.time>parameters.t0+parameters.timespan)) continue;

      // Call the photon detection routine that generates the right charge
      // and stores it in the detector pixels.
      // Before generating and adding the charge to the detector the routine also 
      // checks, whether the integration time is exceeded and performs the readout 
      // in that case. 
      status=addImpact2FramestoreDetector(&detector, &impact);
      if(EXIT_SUCCESS!=status) break;

    } // END of scanning the impact list.

    // Perform a final readout of the FramestoreDetector.
    // Otherwise the last stored charges may be lost.
    status=readoutFramestoreDetector(&detector);
    if(EXIT_SUCCESS!=status) break;

  } while(0); // END of the error handling loop.

  // --- END of Detection process ---


  // --- Cleaning up ---
  headas_chat(5, "\ncleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  if (NULL!=impactlistfile.fptr) fits_close_file(impactlistfile.fptr, &status);

  // Release memory of detector.
  status+=cleanupFramestoreDetector(&detector);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully\n\n");
  return(status);
}




////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int getpar(struct Parameters* parameters)
{
  int status=EXIT_SUCCESS; // error status

  // Get the name of the impact list file (FITS file)
  if ((status = PILGetFname("impactlist_filename", parameters->impactlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the impact list file!\n", status);
  }

  // Get the integration time of the FRAMESTORE CCD.
  else if ((status = PILGetReal("integration_time", &parameters->integration_time))) {
    HD_ERROR_THROW("Error reading the integration time!\n", status);
  }

  // Detector width [pixel]
  else if ((status = PILGetInt("det_width", &parameters->width))) {
    HD_ERROR_THROW("Error reading the width of the detector!\n", status);
  }

  // [m]
  else if ((status = PILGetReal("det_pixelwidth", &parameters->pixelwidth))) {
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

  // Set the event list template file for eROSITA:
  strcat(parameters->eventlist_template, "/erosita.eventlist.tpl");

  return(status);
}


