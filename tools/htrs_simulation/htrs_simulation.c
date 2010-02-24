#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include "htrs_simulation.h"


////////////////////////////////////
/** Main procedure. */
int htrs_simulation_main() {

  // Containing all programm parameters read by PIL
  struct Parameters parameters; 

  // HTRSDetector data structure.
  // This is the main object for the simulation of the HTRS. It contains all
  // important parameters of the detector and is used to create an event file
  // from the incident photon events.
  // Before the first usage it has to be initialized.
  HTRSDetector detector;

  ImpactListFile impactlistfile;

  int status=EXIT_SUCCESS; // Error status.


  // Register HEATOOL:
  set_toolname("htrs_simulation");
  set_toolversion("0.01");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    // Read parameters using PIL library:
    if ((status=getpar(&parameters))) break;
    
    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(SIXT_HD_RANDOM_SEED);


    // Open the impact list FITS file.
    status = openImpactListFile(&impactlistfile, parameters.impactlist_filename, 
				READONLY);
    if (EXIT_SUCCESS!=status) break;


    // Detector settings.
    // Store the settings for the HTRSDetector in the corresponding data structure
    // and call the initialization routine in order to allocate memory etc.
#ifdef HTRS_HEXPIXELS
    struct HTRSDetectorParameters hdparameters = {
      .pixels = { .npixels = 37,
		  .pixelwidth = parameters.pixelwidth },
      .generic = { .ccsigma = parameters.ccsigma, 
		   .pha_threshold = parameters.pha_threshold,
		   .energy_threshold = parameters.energy_threshold,
		   .rmf_filename = parameters.rmf_filename /* String address!! */ },
      .dead_time          = parameters.dead_time,
      .eventlist_filename = parameters.eventlist_filename /* String address!! */,
      .eventlist_template = parameters.eventlist_template /* String address!! */
    };
#endif
#ifdef HTRS_ARCPIXELS
    // Configuration with 31 pixels optimized for uniform photon
    // distribution among the pixels (for photons at 1 keV).
    int npixels[4] = { 1, 6, 12, 12 };
    double radii[4] = { 2.64e-3, 5.5e-3, 8.82e-3, 14.3e-3 };
    double offset_angles[4] = { 0., 0., 0., 0. };
    struct HTRSDetectorParameters hdparameters = {
      .pixels = { .nrings = 4,
		  .npixels = npixels,
		  .radius = radii,
		  .offset_angle = offset_angles,
		  .mask_spoke_width = parameters.mask_spoke_width },
      .generic = { .ccsigma = parameters.ccsigma, 
		   .pha_threshold = parameters.pha_threshold,
		   .energy_threshold = parameters.energy_threshold,
		   .rmf_filename = parameters.rmf_filename /* String address!! */ },
      .dead_time          = parameters.dead_time,
      .eventlist_filename = parameters.eventlist_filename /* String address!! */,
      .eventlist_template = parameters.eventlist_template /* String address!! */
    };
#endif
    if(EXIT_SUCCESS!=(status=initHTRSDetector(&detector, &hdparameters))) break;
    // END of DETECTOR CONFIGURATION SETUP    

    // --- END of Initialization ---


    // --- Beginning of Detection Process ---

    headas_chat(5, "start detection process ...\n");
    Impact impact; // Buffer to store the impacts read from the FITS file.

    // Loop over all impacts in the FITS file.
    while ((EXIT_SUCCESS==status)&&(0==ImpactListFile_EOF(&impactlistfile))) {

      status=getNextImpactListFileRow(&impactlistfile, &impact);
      if (EXIT_SUCCESS!=status) break;

      // Check whether the event lies in the specified time interval:
      if ((impact.time<parameters.t0)||(impact.time>parameters.t0+parameters.timespan)) continue;

      // Call the photon detection routine that generates the right charge
      // and stores it in the detector pixels.
      // Before generating and adding the charge to the detector the routine also 
      // checks, whether the integration time is exceeded and performs the readout 
      // in that case. 
      status=addImpact2HTRSDetector(&detector, &impact);
      if (EXIT_SUCCESS!=status) break;

    } // END of scanning the impact list.
    if (EXIT_SUCCESS!=status) break;
    
  } while(0); // END of the error handling loop.

  // --- END of Detection process ---


  // --- Cleaning up ---
  headas_chat(5, "\ncleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Close the FITS files.
  status += closeImpactListFile(&impactlistfile);

  // Release memory of detector.
  status+=cleanupHTRSDetector(&detector);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully\n\n");
  return(status);
}



////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
static int getpar(struct Parameters* parameters)
{
  int status=EXIT_SUCCESS; // Error status.

  // Get the name of the impact list file (FITS file).
  if ((status = PILGetFname("impactlist_filename", parameters->impactlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the impact list file!\n", status);
  }

  // Get the dead time for a detector pixel.
  else if ((status = PILGetReal("dead_time", &parameters->dead_time))) {
    HD_ERROR_THROW("Error reading the dead time for the detector pixels!\n", status);
  }

#ifdef HTRS_HEXPIXELS
  // [m]
  else if ((status = PILGetReal("pixelwidth", &parameters->pixelwidth))) {
    HD_ERROR_THROW("Error reading the width of the detector pixels!\n", status);
  }
#endif

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

#ifdef HTRS_ARCPIXELS
  // Determine the width of the spokes of the HTRS mask.
  else if ((status = PILGetReal("mask_spoke_width", &parameters->mask_spoke_width))) {
    HD_ERROR_THROW("Error reading the spoke width of the HTRS mask!\n", status);
  }
#endif  

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
  if (EXIT_SUCCESS!=status) return(status);
  // Set the event list template file:
  strcat(parameters->eventlist_template, "/htrs.eventlist.tpl");

  return(status);
}


