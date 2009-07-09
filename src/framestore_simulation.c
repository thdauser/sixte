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
  Detector* detector=NULL;

  struct ImpactlistFile impactlistfile;
  struct Eventlist_File* eventlist_file=NULL;

  int status=EXIT_SUCCESS;   // error status


  // Register HEATOOL:
  set_toolname("framestore_simulation");
  set_toolversion("0.01");


  do {   // beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // Read parameters using PIL library:
    if ((status=getpar(&parameters))) break;
    
    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);

    // Set initial DETECTOR CONFIGURATION.
    detector = get_Detector(&status);
    if (EXIT_SUCCESS!=status) break;

    // General detector settings.
    struct DetectorParameters detectorparameters = {
      .width = parameters.width,
      .pixelwidth = parameters.pixelwidth,
      .ccsigma          = parameters.ccsigma,

      .pha_threshold = parameters.pha_threshold,
      .energy_threshold = parameters.energy_threshold,

      .t0               = parameters.t0
    };
    strcpy(detectorparameters.rmf_filename, parameters.rmf_filename);

    // Framestore-specific settings
    struct FramestoreParameters framestoreparameters = {
      .integration_time = parameters.integration_time,
    };
    init_FramestoreDetector(detector, detectorparameters, framestoreparameters);    
    
    // END of DETECTOR CONFIGURATION SETUP


    // Open the impact list FITS file.
    if(EXIT_SUCCESS!=(status=impactlist_openFile(&impactlistfile, parameters.impactlist_filename, 
						 READONLY))) break;
    strcpy(parameters.attitude_filename, impactlistfile.attitude_filename);
    


    // Delete old EVENT LIST file:
    remove(parameters.eventlist_filename);
    // Create new event list FITS file.
    headas_chat(5, "create FITS file '%s' according to template '%s' ...\n", 
		parameters.eventlist_filename,
		parameters.eventlist_template_filename);
    // Create the new event list file according to the selected template:
    fitsfile* ef_fptr=NULL;
    char buffer[FILENAME_LENGTH];
    sprintf(buffer, "%s(%s/%s)", parameters.eventlist_filename, 
	    parameters.templatedir, parameters.eventlist_template_filename);
    if (fits_create_file(&ef_fptr, buffer, &status)) break;
    if (fits_close_file(ef_fptr, &status)) break;
    
    // Open the newly created event list FITS file for output:
    headas_chat(5, "open FITS file '%s' for output of event list ...\n",
		parameters.eventlist_filename);
    eventlist_file = open_EventlistFile(parameters.eventlist_filename, READWRITE, &status);
    // Add important additional HEADER keywords to the event list.
    if (fits_write_key(eventlist_file->fptr, TSTRING, "ATTITUDE", 
		       parameters.attitude_filename,
		       "name of the attitude FITS file", &status)) break;

    // Set the time-keyword in the Event List Header.
    // See also: Stevens, "Advanced Programming in the UNIX environment", p. 155 ff.
    time_t current_time;
    if (0 != time(&current_time)) {
      struct tm* current_time_utc = gmtime(&current_time);
      if (NULL != current_time_utc) {
	char current_time_str[MAXMSG];
	if (strftime(current_time_str, MAXMSG, "%Y-%m-%dT%H:%M:%S", current_time_utc) > 0) {
	  // Return value should be == 19 !
	  if (fits_update_key(eventlist_file->fptr, TSTRING, "DATE-OBS", current_time_str, 
			     "Start Time (UTC) of exposure", &status)) break;
	}
      }
    } // END of writing time information to Event File FITS header.
    

    // --- END of Initialization ---


    // --- Beginning of Detection Process ---

    headas_chat(5, "start detection process ...\n");
    struct Impact impact, next_impact;
    double next_background_event_time=0.;
    int reached_end_of_impactlist=0;

    // Read the first row of the impact list.
    if(EXIT_SUCCESS!=(status=impactlist_getNextRow(&impactlistfile, &next_impact))) break;

    while (EXIT_SUCCESS==status) {
      // TODO: Break the loop, when interval time+timespan is exceeded.

      if ((parameters.background_rate>0.) && (next_background_event_time<next_impact.time)) {
	// The current event is a background event:
	impact.time = next_background_event_time;
	impact.energy = 1.; // TODO
	impact.position.x = 2*(get_random_number()-0.5)*(detector->offset*detector->pixelwidth);
	impact.position.y = 2*(get_random_number()-0.5)*(detector->offset*detector->pixelwidth);
	// TODO: prevent PSF check for these events !!
	
	// Determine the time of the NEXT background event:
	next_background_event_time += rndexp(1./(double)parameters.background_rate);

      } else {
	// The current event is obtained from the impact list.
	impact = next_impact;

	if (0==reached_end_of_impactlist) {
	  // Read in the next row from the impact list:
	  if(EXIT_SUCCESS!=(status=impactlist_getNextRow(&impactlistfile, &next_impact))) break;
	  // Check if we reached the end of the impact list:
	  if (impactlistfile.row >= impactlistfile.nrows) reached_end_of_impactlist = 1;
	} else {
	  // There are no further rows in the impact list. So we have to stop here.
	  break;
	}
      }

      // Call the detector action routine: this routine checks, whether the 
      // integration time is exceeded and performs the readout in this case. 
      // Otherwise it will simply do nothing.
      if (detector->readout != NULL) { // HTRS and TES do not have this routine.
	detector->readout(detector, impact.time, eventlist_file, &status);
	if (status != EXIT_SUCCESS) break;
      }

      // Check whether the event lies in the specified time interval:
      if ((impact.time > parameters.t0) && (impact.time < parameters.t0+parameters.timespan)) {

	// Measure the photon in the detector pixels, i.e., create the 
	// corresponding charges there.

	// Determine a detector channel (PHA channel) according to RMF.
	// The channel is obtained from the RMF using the corresponding
	// HEAdas routine which is based on drawing a random number.
	long channel;
	ReturnChannel(detector->rmf, impact.energy, 1, &channel);

	// Check if the photon is really measured. If the
	// PHA channel returned by the HEAdas RMF function is '-1', 
	// the photon is not detected.
	if (channel==-1) {
	  continue;  // -> Continue with the next photon in the list.
	}
	assert(channel>=0);

	// Get the corresponding created charge.
	// NOTE: In this simulation the charge is represented by the nominal
	// photon energy which corresponds to the PHA channel according to the
	// EBOUNDS table.
	float charge = get_energy(channel, detector);
	
	if(charge > 0.) {
	  int x[4], y[4];
	  double fraction[4];
	  
	  // Determine the affected detector pixels.
	  int npixels = get_pixel_square(detector, impact.position, x, y, fraction);

	  // Add the charge created by the photon to the affected detector pixels.
	  int count;
	  for (count=0; count<npixels; count++) {
	    if (x[count] != INVALID_PIXEL) {
	      detector->pixel[x[count]][y[count]].charge += 
		charge * fraction[count] * 
		// |      |-> charge fraction due to split events
		// |-> charge created by incident photon
		detector_active(x[count], y[count], detector, impact.time);
	      // |-> "1" if pixel can measure charge, "0" else
	    }
	  }
	  
	} // END if(charge>0.)
      } // END 'time' within specified time interval
    } // END of scanning the impact list.

  } while(0);  // END of the error handling loop.

  // --- END of Detection process ---


  // --- Cleaning up ---
  headas_chat(5, "\ncleaning up ...\n");

  // Release HEADAS random number generator
  HDmtFree();

  if (NULL!=impactlistfile.fptr) fits_close_file(impactlistfile.fptr, &status);

  if (NULL!=eventlist_file) {
    if (eventlist_file->fptr) fits_close_file(eventlist_file->fptr, &status);
  }

  // Release memory of detector:
  if (detector != NULL) {
    if (detector->pixel != NULL) {
      int count;
      for (count = 0; count < detector->width; count++) {
	if (detector->pixel[count]) {
	  free(detector->pixel[count]);
	}
      }
      free(detector->pixel);
    }
  }

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
      strcpy(parameters->templatedir, buffer);
    } else {
      if ((status = PILGetFname("fits_templates", parameters->templatedir))) {
	HD_ERROR_THROW("Error reading the path of the FITS templates!\n", status);
      }
    }
  }

  // Set the event list template file for eROSITA:
  strcpy(parameters->eventlist_template_filename, "erosita.eventlist.tpl");

  return(status);
}


