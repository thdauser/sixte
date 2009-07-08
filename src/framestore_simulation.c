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

  fitsfile *impactlist_fptr=NULL; 
  long impactlist_nrows;

  // Detector data structure (containing the pixel array, its width, ...)
  Detector* detector=NULL;

  struct Eventlist_File* eventlist_file=NULL;

  char msg[MAXMSG];          // error output buffer
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

    // GENERAL SETTINGS
    detector->type = FRAMESTORE;
    detector->width = parameters.width;
    detector->offset = detector->width/2;
    detector->pixelwidth = parameters.pixelwidth;

    // Event thresholds:
    detector->pha_threshold = parameters.pha_threshold;
    detector->energy_threshold = parameters.energy_threshold;

    // DETECTOR SPECIFIC SETTINGS
    struct FramestoreParameters framestoreparameters = {
      .integration_time = parameters.integration_time,
      .ccsigma          = parameters.ccsigma
    };
    init_FramestoreDetector(detector, framestoreparameters);    
    
    // Get the memory for the detector pixels
    if (get_DetectorPixels(detector, &status)) break;

    // Set the first readout time such that the first readout is performed 
    // immediately at the beginning of the simulation.
    detector->readout_time = parameters.t0; // TODO move to init_FramestoreDetector() routine.
      
     
    // Read the detector RMF and EBOUNDS from the specified file and 
    // assign them to the Detector data structure.
    if ((status=detector_assign_rsp(detector, parameters.rmf_filename)) 
	!= EXIT_SUCCESS) break;


    // END of DETECTOR CONFIGURATION SETUP


    // Open impact list FITS file:
    if (fits_open_table(&impactlist_fptr, parameters.impactlist_filename, 
			READONLY, &status)) break;
    // Determine the number of rows in the impact list:
    if (fits_get_num_rows(impactlist_fptr, &impactlist_nrows, &status)) break;
    if (0 >= impactlist_nrows) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: impact list in file '%s' is empty!\n", 
	      parameters.impactlist_filename);
      HD_ERROR_THROW(msg, status);
      break;
    }

    // Read HEADER keywords.
    char comment[MAXMSG]; // buffer
    if (fits_read_key(impactlist_fptr, TSTRING, "ATTITUDE", 
		      &parameters.attitude_filename, 
		      comment, &status)) break;


    // Delete old event list file:
    remove(parameters.eventlist_filename);
    // Create new event list FITS file.
    headas_chat(5, "create FITS file '%s' according to template '%s' ...\n", 
		parameters.eventlist_filename,
		parameters.eventlist_template_filename);
    // Create the new event list file according to the selected template:
    fitsfile* ef_fptr=NULL;
    char buffer[FILENAME_LENGTH];
    sprintf(buffer, "%s(%s)", parameters.eventlist_filename, 
	    parameters.eventlist_template_filename);
    if (fits_create_file(&ef_fptr, buffer, &status)) break;
    if (fits_close_file(ef_fptr, &status)) break;
    
    // Open the newly created event list FITS file for output:
    headas_chat(5, "open FITS file '%s' for output of event list ...\n",
		parameters.eventlist_filename);
    eventlist_file = open_EventlistFile(parameters.eventlist_filename, READWRITE, &status);
    // Create a new FITS file and a table for the event list:


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

    int anynul;
    double time, next_real_impact_time=0., next_background_event_time=0.;
    float energy;
    struct Point2d position;
    long impactlist_row=0;
    while (status==EXIT_SUCCESS) {
      // TODO: Break the loop, when interval time+timespan is exceeded.

      if ((parameters.background_rate > 0.) && 
	  (next_background_event_time < next_real_impact_time)) {
	// The current event is a background event:
	time = next_background_event_time;
	energy = 1.; // TODO
	position.x = 2*(get_random_number()-0.5) * (detector->offset*detector->pixelwidth);
	position.y = 2*(get_random_number()-0.5) * (detector->offset*detector->pixelwidth);
	// TODO: prevent PSF check for these events !!
	
	// Determine the time of the NEXT background event:
	next_background_event_time += rndexp(1./(double)parameters.background_rate);

      } else {
	// The current event is obtained from the impact list:
	// Read an entry from the impact list:
	anynul = 0;
	time = next_real_impact_time;
	energy = 0.;
	position.x = 0.;
	position.y = 0.;
	fits_read_col(impactlist_fptr, TFLOAT, 2, impactlist_row+1, 1, 1, 
		      &energy, &energy, &anynul, &status);
	fits_read_col(impactlist_fptr, TDOUBLE, 3, impactlist_row+1, 1, 1, 
		      &position.x, &position.x, &anynul, &status);
	fits_read_col(impactlist_fptr, TDOUBLE, 4, impactlist_row+1, 1, 1, 
		      &position.y, &position.y, &anynul, &status);

	// Get the time of the next entry in the impact list:
	impactlist_row++;
	if (impactlist_row>=impactlist_nrows) break; // Reached end of impact list.
	next_real_impact_time = 0.;
	fits_read_col(impactlist_fptr, TDOUBLE, 1, impactlist_row+1, 1, 1, 
		      &next_real_impact_time, &next_real_impact_time, &anynul, &status);
      }

      // Call the detector action routine: this routine checks, whether the 
      // integration time is exceeded and performs the readout in this case. 
      // Otherwise it will simply do nothing.
      if (detector->readout != NULL) { // HTRS and TES do not have this routine.
	detector->readout(detector, time, eventlist_file, &status);
	if (status != EXIT_SUCCESS) break;
      }

      // Check whether the event lies in the specified time interval:
      if ((time > parameters.t0) && (time < parameters.t0+parameters.timespan)) {

	// Measure the photon in the detector pixels, i.e., create the 
	// corresponding charges there.

	// Determine a detector channel (PHA channel) according to RMF.
	// The channel is obtained from the RMF using the corresponding
	// HEAdas routine which is based on drawing a random number.
	long channel;
	ReturnChannel(detector->rmf, energy, 1, &channel);

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
      
	  if ((detector->type == FRAMESTORE) || (detector->type == WFI)) {
	    // Determine the affected detector pixels.
	    int npixels = get_pixel_square(detector, position, x, y, fraction);

	    // Add the charge created by the photon to the affected detector pixels.
	    int count;
	    for (count=0; count<npixels; count++) {
	      if (x[count] != INVALID_PIXEL) {
		detector->pixel[x[count]][y[count]].charge += 
		  charge * fraction[count] * 
		  // |      |-> charge fraction due to split events
		  // |-> charge created by incident photon
		  detector_active(x[count], y[count], detector, time);
		// |-> "1" if pixel can measure charge, "0" else
	      }
	    }
	    
	  } else if (detector->type == HTRS) {
	    struct Point2d position2;
	    position2.x = position.y;
	    position2.y = position.x; // TODO !!!!!!!!!!!!
	    int npixels = htrs_get_pixel(detector, position2, x, y, fraction);
	    
	    struct Event event;
	    int count;
	    for (count=0; count<npixels; count++) {
	      if (x[count] != INVALID_PIXEL) {
		// Check if the affected detector pixel is active:
		if (htrs_detector_active(x[count], y[count], detector, time)) {
		  
		  // Save the time of the photon arrival in this pixel
		  detector->pixel[x[count]][y[count]].arrival = time;
		
		  // Store the photon charge and the new arrival time:
		  event.pha = get_channel(charge*fraction[count], detector);
		  event.time = time;                // TODO: drift time
		  event.xi = detector->htrs_icoordinates2pixel[x[count]][y[count]]+1;
		  event.yi = 0;  // human readable HTRS pixels numbers start at 1 <-|
		  event.frame = 0;
		  event.ra = NAN;
		  event.dec = NAN;
		  event.sky_xi = 0;
		  event.sky_yi = 0;
		  
		  // Add the event to the FITS event list.
		  // Check lower PHA threshold:
		  if ((event.pha >= detector->pha_threshold)&&
		      (charge*fraction[count] >= detector->energy_threshold)){ 
		    // There is an event in this pixel, so insert it into eventlist:
		    add_eventlist_row(eventlist_file, event, &status);
		  }
		} // END htrs_detector_active(...)
	      } // END x[count] != INVALID_PIXEL
	    } // END of loop over all split partners.
	
	  } else if (detector->type == TES) {
	    get_pixel_square(detector, position, x, y, fraction);
	    
	    if (x[0] != INVALID_PIXEL) {
	      struct Event event;
	    
	      // Store the photon charge and the new arrival time:
	      event.pha = get_channel(charge, detector);  // TODO: RMF
	      event.time = time;
	      event.xi = x[0];
	      event.yi = y[0];
	      event.frame = detector->frame;
	      event.ra = NAN;
	      event.dec = NAN;
	      event.sky_xi = 0;
	      event.sky_yi = 0;
	      
	      // Add the event to the FITS event list.
	      if ((event.pha>=detector->pha_threshold)&&
		  (charge>=detector->energy_threshold)){ // Check lower PHA threshold
		// There is an event in this pixel, so insert it into eventlist:
		add_eventlist_row(eventlist_file, event, &status);
	      }
	    } // END x[0] != INVALID_PIXEL
	  } // END detector->type == TES
	} // END if(charge>0.)
      } // END 'time' within specified time interval
    } // END of scanning the impact list.

  } while(0);  // END of the error handling loop.

  // --- END of Detection process ---


  // --- Cleaning up ---
  headas_chat(5, "\ncleaning up ...\n");

  // Release HEADAS random number generator
  HDmtFree();

  if (impactlist_fptr) fits_close_file(impactlist_fptr, &status);
  if (NULL!=eventlist_file) {
    if (eventlist_file->fptr) fits_close_file(eventlist_file->fptr, &status);
  }

  // Release memory of detector:
  int count;
  if (detector->pixel) {
    for (count = 0; count < detector->width; count++) {
      if (detector->pixel[count]) {
	free(detector->pixel[count]);
      }
    }
    free(detector->pixel);
  }

  if (detector->type == HTRS) {
    htrs_free_Detector(detector);
  }

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully\n\n");

  return(status);
}




////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int getpar(struct Parameters* parameters)
{
  char msg[MAXMSG];             // error output buffer
  int status=EXIT_SUCCESS;      // error status

  // Get the name of the impact list file (FITS file)
  if ((status = PILGetFname("impactlist_filename", parameters->impactlist_filename))) {
    sprintf(msg, "Error reading the name of the impact list file!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the integration time of the FRAMESTORE CCD.
  else if ((status = PILGetReal("integration_time", &parameters->integration_time))) {
    sprintf(msg, "Error reading the integration time!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Detector width [pixel]
  else if ((status = PILGetInt("det_width", &parameters->width))) {
    sprintf(msg, "Error reading the width of the detector!\n");
    HD_ERROR_THROW(msg,status);
  }

  // [m]
  else if ((status = PILGetReal("det_pixelwidth", &parameters->pixelwidth))) {
    sprintf(msg, "Error reading the width of the detector pixels!\n");
    HD_ERROR_THROW(msg,status);
  }

  // [m]
  else if ((status = PILGetReal("ccsigma", &parameters->ccsigma))) {
    sprintf(msg, "Error reading the charge cloud sigma!\n");
    HD_ERROR_THROW(msg,status);
  }
  if (status) return(status);

  
  // Read the detector thresholds (either integer PHA or float energy):
  int pha_threshold;
  if ((status = PILGetInt("pha_threshold", &pha_threshold))) {
    sprintf(msg, "Error: could not determine detector PHA threshold!\n");
    HD_ERROR_THROW(msg,status);
    return(status);
  } else {
    parameters->pha_threshold = (long)pha_threshold;
  }
  if (parameters->pha_threshold==-1) {
    if ((status = PILGetReal4("energy_threshold", &parameters->energy_threshold))) {
      sprintf(msg, "Error: could not determine detector energy threshold!\n");
      HD_ERROR_THROW(msg,status);
      return(status);
    }
  } else {
    parameters->energy_threshold=0.;
  }

  // Get the name of the detector redistribution file (FITS file)
  if ((status = PILGetFname("rmf_filename", parameters->rmf_filename))) {
    sprintf(msg, "Error reading the name of the detector" 
	    "redistribution matrix file (RMF)!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the background count rate
  else if ((status = PILGetReal4("background_rate", &parameters->background_rate))) {
    sprintf(msg, "Error: could not determine the detector background rate!\n");
    HD_ERROR_THROW(msg,status);
    return(status);
  }

  // Get the name of the output event list (FITS file)
  else if ((status = PILGetFname("eventlist_filename", parameters->eventlist_filename))) {
    sprintf(msg, "Error reading the name of the event list file!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the name of the output event list TEMPLATE (ASCII file)
  else if ((status = PILGetFname("eventlist_template_filename", 
				 parameters->eventlist_template_filename))) {
    HD_ERROR_THROW("Error reading the name of the event list file!\n", status);
  }

  // Get the start time of the simulation
  else if ((status = PILGetReal("t0", &parameters->t0))) {
    sprintf(msg, "Error reading the 't0' parameter!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the timespan for the simulation
  else if ((status = PILGetReal("timespan", &parameters->timespan))) {
    sprintf(msg, "Error reading the 'timespan' parameter!\n");
    HD_ERROR_THROW(msg,status);
  }

  return(status);
}


