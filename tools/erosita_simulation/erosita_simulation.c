#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include "erosita_simulation.h"


////////////////////////////////////
/** Main procedure. */
int erosita_simulation_main() {
  struct Parameters parameters; // Containing all programm parameters read by PIL

  // Detector data structure (containing the pixel array, its width, ...).
  eROSITADetectors detector;

  // Data structure to model the detector background.
  UniformDetectorBackground background;

  ImpactListFile impactlistfile;
  // WCS reference values for the position of the orginial input data.
  // These data are needed for the eROSITA image reconstruction algorithm
  // in order to determine the right WCS header keywords for, e.g., cluster images.
  double refxcrvl, refycrvl; 
  // Filename of the attitude file.
  char attitude_filename[MAXMSG];

  int status=EXIT_SUCCESS; // Error status.


  // Register HEATOOL:
  set_toolname("erosita_simulation");
  set_toolversion("0.01");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    // Read parameters using PIL library:
    if ((status=getpar(&parameters))) break;


    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);


    // DETECTOR setup.
    struct eROSITADetectorsParameters fdparameters = {
      .pixels = 
      { .xwidth = parameters.width,
	.ywidth = parameters.width,
	.xpixelwidth = parameters.pixelwidth,
	.ypixelwidth = parameters.pixelwidth 
      },
      .generic = 
      { .ccsigma          = parameters.ccsigma, 
	.pha_threshold    = parameters.pha_threshold,
	.energy_threshold = parameters.energy_threshold,
	.rmf_filename     = parameters.rmf_filename /* String address!! */ 
      },
      .integration_time   = parameters.integration_time,
      .t0                 = parameters.t0,
      .eventlist_filename = parameters.eventlist_filename /* String address!! */,
      .eventlist_template = parameters.eventlist_template
    };
    status=initeROSITADetectors(&detector, &fdparameters);
    if(EXIT_SUCCESS!=status) break;
    // END of DETECTOR CONFIGURATION SETUP    

    // Open the FITS file with the input impact list:
    status = openImpactListFile(&impactlistfile, parameters.impactlist_filename, 
				READONLY);
    if (EXIT_SUCCESS!=status) break;
    // Determine the WCS keywords.
    char comment[MAXMSG]; // buffer
    if (fits_read_key(impactlistfile.fptr, TDOUBLE, "REFXCRVL", &refxcrvl, 
		      comment, &status)) break;    
    if (fits_read_key(impactlistfile.fptr, TDOUBLE, "REFYCRVL", &refycrvl, 
		      comment, &status)) break;    
    // Determine the name of the attitude file from the FITS header.
    if (fits_read_key(impactlistfile.fptr, TSTRING, "ATTITUDE", 
		      attitude_filename, comment, &status)) break;

    // Add important additional HEADER keywords to the event list.
    char keyword[MAXMSG];
    sprintf(keyword, "TCRVL%d", detector.eventlist.cskyx);
    if (fits_update_key(detector.eventlist.generic.fptr, TDOUBLE, keyword, &refxcrvl, 
			"", &status)) break;
    sprintf(keyword, "TCRVL%d", detector.eventlist.cskyy);
    if (fits_update_key(detector.eventlist.generic.fptr, TDOUBLE, keyword, &refycrvl, 
			"", &status)) break;
    if (fits_update_key(detector.eventlist.generic.fptr, TSTRING, "ATTITUDE", 
			attitude_filename, "name of the attitude FITS file", &status)) break;


    // Setup detector BACKGROUND data structure.
    struct UniformDetectorBackgroundParameters bkgdparameters = {
      .rate              = parameters.background_rate,
      .spectrum_filename = parameters.background_filename /* String address!! */
    };
    if(EXIT_SUCCESS!=(status=initUniformDetectorBackground(&background, &bkgdparameters))) break;
    // END of BACKGROUND CONFIGURATION SETUP

    // --- END of Initialization ---



    // --- Beginning of Detection Process ---

    headas_chat(5, "start detection process ...\n");
    Impact impact;
    double former_impact_time=-1000.; // Time of the last impact.

    // Loop over all impacts in the FITS file.
    while ((EXIT_SUCCESS==status)&&(0==ImpactListFile_EOF(&impactlistfile))) {

      status=getNextImpactListFileRow(&impactlistfile, &impact);
      if(EXIT_SUCCESS!=status) break;

      // Check whether the event lies in the specified time interval:
      if ((impact.time<parameters.t0)||(impact.time>parameters.t0+parameters.timespan)) 
	continue;

      if (background.rate > 0.) { // Insert Background events.
	// If the difference between to subsequent cosmic photon impacts
	// is greater than 60s, interrupt the background generation in order
	// to avoid numerical costs during periods in which the telescope is not
	// pointing to any interesting part in the sky. These times will be skipped
	// anyway.
	if (impact.time-former_impact_time > 60.) {
	  // Fill up the time in the 60s after the last photon impact with background events.
	  while (background.nextImpact.time-former_impact_time<60.) {
	    // Add the background event to the CCD array.
	    status=addImpact2eROSITADetectors(&detector, &background.nextImpact);
	    if(EXIT_SUCCESS!=status) break;
	    // Create a new background event.
	    createUniformDetectorBackgroundImpact(&background, &detector.pixels[0],
						  detector.generic.rmf);
	  }
	  if(EXIT_SUCCESS!=status) break;
	  // Jump to the next interval where cosmic photons arrive at the detector.
	  // Leave out the time in between in order to reduce numerical costs.
	  Impact nextImpact = { .time = impact.time - 60.,
				.energy = 0.,
				.position = { .x=0., .y=0. } };
	  background.nextImpact = nextImpact;
	}
      
	// Insert uniformly distributed detector background events in the interval between
	// the last background event and the current photon arrival time.
	while (background.nextImpact.time < impact.time) {
	  // Add the background event to the CCD array.
	  status=addImpact2eROSITADetectors(&detector, &background.nextImpact);
	  if(EXIT_SUCCESS!=status) break;
	  // Create a new background event.
	  createUniformDetectorBackgroundImpact(&background, &detector.pixels[0], 
						detector.generic.rmf);
	}
	if(EXIT_SUCCESS!=status) break;
      }	// END of inserting background events.

      // Call the photon detection routine that generates the right charge
      // and stores it in the detector pixels.
      // Before generating and adding the charge to the detector the routine also 
      // checks, whether the integration time is exceeded and performs the readout 
      // in that case. 
      status=addImpact2eROSITADetectors(&detector, &impact);
      if(EXIT_SUCCESS!=status) break;

      // Save the time of the current impact (necessary for background generation).
      former_impact_time = impact.time;

    } // END of scanning the impact list.

    // Perform a final readout of the eROSITADetectors.
    // Otherwise the last stored charges may be lost.
    status=readouteROSITADetectors(&detector);
    if(EXIT_SUCCESS!=status) break;

  } while(0); // END of the error handling loop.

  // --- END of Detection process ---


  // --- Cleaning up ---
  headas_chat(5, "\ncleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Close the FITS files.
  status += closeImpactListFile(&impactlistfile);

  // Release memory of detector.
  status+=cleanupeROSITADetectors(&detector);
  
  // Release memory of background data structure.
  status+=cleanupUniformDetectorBackground(&background);

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

  // Get the integration time of the eROSITA CCDs.
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

  // Get the filename of the detector background spectrum.
  else if ((status = PILGetFname("background_filename", parameters->background_filename))) {
    HD_ERROR_THROW("Error: could not determine the detector background spectrum!\n", status);
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

