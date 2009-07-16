#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include "wfi_simulation.h"


////////////////////////////////////
/** Main procedure. */
int wfi_simulation_main() {
  struct Parameters parameters; // Containing all programm parameters read by PIL

  // Detector data structure (containing the pixel array, its width, ...)
  WFIDetector detector;

  struct ImpactlistFile impactlistfile;

  int status=EXIT_SUCCESS; // error status


  // Register HEATOOL:
  set_toolname("wfi_simulation");
  set_toolversion("0.01");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    // Read parameters using PIL library:
    if ((status=getpar(&parameters))) break;
    
    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);


    // Open the impact list FITS file.
    if(EXIT_SUCCESS!=(status=impactlist_openFile(&impactlistfile, parameters.impactlist_filename, 
						 READONLY))) break;

    // Detector settings.
    struct WFIDetectorParameters fdparameters = {
      .pixels = { .xwidth = parameters.width,
		  .ywidth = parameters.width,
		  .xpixelwidth = parameters.pixelwidth,
		  .ypixelwidth = parameters.pixelwidth },
      .generic = { .ccsigma = parameters.ccsigma, 
		   .pha_threshold = parameters.pha_threshold,
		   .energy_threshold = parameters.energy_threshold,
		   .rmf_filename = parameters.rmf_filename /* String address!! */ },
      .readout_directions = parameters.readout_directions,
      .line_readout_time  = parameters.line_readout_time,
      .line_clear_time    = parameters.line_clear_time,
      .t0                 = parameters.t0,
      .eventlist_filename = parameters.eventlist_filename /* String address!! */,
      .eventlist_template = parameters.eventlist_template /* String address!! */
    };
    if(EXIT_SUCCESS!=(status=initWFIDetector(&detector, &fdparameters))) break;

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

      // Call the detector readout routine: this routine checks, whether the 
      // integration time is exceeded and performs the readout in this case. 
      // Otherwise it will simply do nothing.
      status=checkReadoutWFIDetector(&detector, impact.time);
      if (EXIT_SUCCESS!=status) break;
	
      // Call the photon detection routine that generates the right charge
      // and stores it in the detector pixels.
      addImpact2WFIDetector(&detector, &impact);

    } // END of scanning the impact list.

    // TODO: Perform one last readout !!
    // Otherwise charges/impacts might be lost.
    status=checkReadoutWFIDetector(&detector, parameters.t0+parameters.timespan);
    if (EXIT_SUCCESS!=status) break;
    
  } while(0); // END of the error handling loop.

  // --- END of Detection process ---


  // --- Cleaning up ---
  headas_chat(5, "\ncleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  if (NULL!=impactlistfile.fptr) fits_close_file(impactlistfile.fptr, &status);

  // Release memory of detector.
  cleanupWFIDetector(&detector);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully\n\n");
  return(status);
}




////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
static int getpar(struct Parameters* parameters)
{
  int status=EXIT_SUCCESS; // error status

  // Get the name of the impact list file (FITS file)
  if ((status = PILGetFname("impactlist_filename", parameters->impactlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the impact list file!\n", status);
  }

  // Number of readout directions.
  else if ((status = PILGetInt("readout_directions", &parameters->readout_directions))) {
    HD_ERROR_THROW("Error reading the number of readout directions!\n", status);
  }

  // Get the readout time for one detector line.
  else if ((status = PILGetReal("line_readout_time", &parameters->line_readout_time))) {
    HD_ERROR_THROW("Error reading the readout time per detector line!\n", status);
  }

  // Get the clear time for one detector line.
  else if ((status = PILGetReal("line_clear_time", &parameters->line_clear_time))) {
    HD_ERROR_THROW("Error reading the clear time per detector line!\n", status);
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

  // Set the event list template file for the different WFI modes:
  char template_filename[MAXMSG];
  if (16==parameters->width) {
    strcpy(template_filename, "wfi.window16.eventlist.tpl");
  } else if (1024==parameters->width) {
    strcpy(template_filename, "wfi.full1024.eventlist.tpl");
  } else {
    status = EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "Error: detector width (%d pixels) is not supported!\n", parameters->width);
    HD_ERROR_THROW(msg, status);
    return(status);
  }
  strcat(strcat(parameters->eventlist_template, "/"), template_filename);

  return(status);
}


