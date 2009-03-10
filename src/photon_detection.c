#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include "photon_detection.h"



////////////////////////////////////
// Main procedure.
int photon_detection_main() {
  double t0;        // starting time of the simulation [s]
  double timespan;  // time span of the simulation [s]

  char impactlist_filename[FILENAME_LENGTH]; // input: impact list
  fitsfile *impactlist_fptr=NULL; 
  long impactlist_nrows;
  char rmf_filename[FILENAME_LENGTH];        // input: RMF

  // Detector data structure (containing the pixel array, its width, ...)
  Detector* detector;   

  struct Eventlist_File eventlist_file;

  char msg[MAXMSG];          // error output buffer
  int status=EXIT_SUCCESS;   // error status


  // Register HEATOOL:
  set_toolname("photon_detection");
  set_toolversion("0.01");

  do {   // beginning of the error handling loop (will at most be run once)

    // --- Initialization ---

    // Read parameters using PIL library:
    detector = get_Detector(&status);
    if (status!=EXIT_SUCCESS) break;

    if ((status=photon_detection_getpar(impactlist_filename, rmf_filename,
					eventlist_file.filename, &t0, &timespan,
					detector))) break;


    // Initialize HEADAS random number generator and GSL generator for 
    // Gaussian distribution.
    HDmtInit(1);


    // Set initial DETECTOR CONFIGURATION.

    // GENERAL SETTINGS

    // Get the memory for the detector pixels
    if (get_DetectorPixels(detector, &status)) break;

    // Calculate initial parameter values from the PIL parameters:
    detector->offset = detector->width/2;

    // Size of the charge cloud [real pixel]
    detector->ccsize = 3.*detector->ccsigma;

    // Set the current detector frame to its initial value:
    detector->frame = -1;
    
    // DETECTOR SPECIFIC SETTINGS
    if (detector->type == FRAMESTORE) {
      headas_chat(5, "--> FRAMESTORE <--\n");

      // Set the first readout time such that the first readout is performed 
      // immediately at the beginning of the simulation (FRAMESTORE).
      detector->readout_time = t0;
      detector->frame = 0;
      
      detector->action = framestore_detector_action;

    } else if (detector->type == DEPFET) {
      headas_chat(5, "--> DEPFET <--\n");
      
      // Set the first readout time such that the first readout is performed 
      // immediately at the beginning of the simulation (DEPFET).
      detector->readout_time = t0 - detector->dead_time; 
      // The readout process starts at the center of the WFI detector, 
      // but for that purpose the current line has to be set to 0, so that the
      // first readout is performed in the middle of the detector array.
      detector->readout_line = 0;

      detector->readout_line = detector->width;   // REMOVE !!!
      detector->action = depfet_detector_action2; // !!!
      
      headas_chat(5, "dead time: %lf\n", detector->dead_time);
      headas_chat(5, "clear time: %lf\n", detector->clear_time);
      headas_chat(5, "readout time: %lf\n", detector->readout_time);
      
    } else if (detector->type == TES) {
      headas_chat(5, "--> TES Microcalorimeter Array <--\n");
    
      detector->action = NULL; // tes_detector_action;

    } else if (detector->type == HTRS) {
      headas_chat(5, "--> HTRS <--\n");

      detector->action = NULL; // htrs_detector_action;

    } else {

      status=EXIT_FAILURE;
      sprintf(msg, "Error: invalid detector type!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }
    if (status != EXIT_SUCCESS) break;
    
    // Consistency check for size of charge cloud:
    if (detector->ccsize > 1.) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: charge cloud size greater than pixel width!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }

    // Get the energy bins of the PHA channels
    if ((status=get_ebounds(&detector->ebounds, &detector->Nchannels, rmf_filename))
	!=EXIT_SUCCESS) break;

    // Get the detector redistribution matrix (RMF)
    if ((status=get_rmf(detector, rmf_filename)) != EXIT_SUCCESS) break;

    // Print some debug information:
    headas_chat(5, "detector pixel width: %lf m\n", detector->pixelwidth);
    headas_chat(5, "charge cloud size: %lf m\n", detector->ccsize);
    headas_chat(5, "number of PHA channels: %d\n", detector->Nchannels);
    headas_chat(5, "PHA threshold: %ld\n", detector->pha_threshold);
    headas_chat(5, "energy threshold: %lf\n\n", detector->energy_threshold);

    // END of DETECTOR CONFIGURATION SETUP


    // Open impact list FITS file:
    if (fits_open_table(&impactlist_fptr, impactlist_filename, 
			READONLY, &status)) break;
    // Determine the number of rows in the impact list:
    if (fits_get_num_rows(impactlist_fptr, &impactlist_nrows, &status)) break;


    // Create event list FITS file:
    headas_chat(5, "create FITS file '%s' for output of event list ...\n", 
		eventlist_file.filename);
    // Delete old event list:
    remove(eventlist_file.filename);
    // Create a new FITS file and a table for the event list:
    if (create_eventlist_file(&eventlist_file, detector, t0, t0+timespan, 
		        // HEADER keywords for event list FITS file:
			// TELESCOP    CCD       INSTRUME
			  "eROSITA",  "pnCCD1", "eROSITA",  &status)) break;

    // --- END of Initialization ---


    // --- Beginning of Detection Process ---

    headas_chat(5, "start detection process ...\n");

    long impactlist_row;
    for(impactlist_row=0; (impactlist_row<impactlist_nrows)&&(status==EXIT_SUCCESS); 
	impactlist_row++) {

      // Read an entry from the impact list:
      int anynul = 0;
      double time = 0.;
      float energy = 0.;
      struct Point2d position;
      position.x = 0.;
      position.y = 0.;
      fits_read_col(impactlist_fptr, TDOUBLE, 1, impactlist_row+1, 1, 1, 
		    &time, &time, &anynul, &status);
      fits_read_col(impactlist_fptr, TFLOAT, 2, impactlist_row+1, 1, 1, 
		    &energy, &energy, &anynul, &status);
      fits_read_col(impactlist_fptr, TDOUBLE, 3, impactlist_row+1, 1, 1, 
		    &position.x, &position.x, &anynul, &status);
      fits_read_col(impactlist_fptr, TDOUBLE, 4, impactlist_row+1, 1, 1, 
		    &position.y, &position.y, &anynul, &status);


      // Call the detector action routine: this routine checks, whether the 
      // integration time is exceeded and performs the readout in this case. 
      // Otherwise it will simply do nothing.
      detector->action(detector, time, &eventlist_file, &status);
      if (status != EXIT_SUCCESS) break;


      // Check whether the event lies in the specified time interval:
      if ((time > t0) && (time < t0+timespan)) {

	// Measure the photon in the detector pixels, i.e., create the 
	// corresponding charges there.

	// Determine a detector channel (PHA channel) according to RMF.
	long channel = detector_rmf(energy, &detector->rmf);
	// Get the corresponding created charge.
	float charge = get_charge(channel, &detector->ebounds);
      
	assert(channel >= 0);   
	assert(channel <= 4096);

	int x[4], y[4];
	double fraction[4];
      
	if ((detector->type == FRAMESTORE) || (detector->type == DEPFET)) {
	  // Determine the affected detector pixels.
	  int npixels = get_pixel_square(detector, position, x, y, fraction);

	  // Add the charge created by the photon to the affected detector pixels.
	  int count;
	  for (count=0; count<npixels; count++) {
	    if (x[count] != INVALID_PIXEL) {
	      detector->pixel[x[count]][y[count]].charge += 
		charge * fraction[count] * 
		// |        |-> charge fraction due to split events
		// |-> charge created by incident photon
		detector_active(x[count], y[count], detector, time);
	      // |-> "1" if pixel can measure charge, "0" else
	    }
	  }

	} else if (detector->type == HTRS) {
	  int npixels = htrs_get_pixel(detector, position, x, y, fraction);

	  struct Event event;
	  int count;
	  for (count=0; count<npixels; count++) {
	    if (x[count] != INVALID_PIXEL) {
	      // Check if the affected detector pixel is active:
	      if (htrs_detector_active(x[count], y[count], detector, time)) {
		
		// Save the time of the photon arrival in this pixel
		detector->pixel[x[count]][y[count]].arrival = time;
		
		// Store the photon charge and the new arrival time:
		event.pha = get_pha(charge*fraction[count], detector);
		event.time = time;                // TODO: drift time
		event.xi = detector->htrs_icoordinates2pixel[x[count]][y[count]]+1;
		event.yi = 0;  // human readable HTRS pixels numbers start at 1 <-|
		event.grade = 0;
		event.frame = 0;

		// Add the event to the FITS event list.
		// Check lower PHA threshold:
		if ((event.pha >= detector->pha_threshold)&&
		    (charge*fraction[count] >= detector->energy_threshold)){ 
		  // There is an event in this pixel, so insert it into eventlist:
		  add_eventlist_row(&eventlist_file, event, &status);
		}
	      } // END htrs_detector_active(...)
	    } // END x[count] != INVALID_PIXEL
	  } // END of loop over all split partners.
	
	} else if (detector->type == TES) {
	  get_pixel_square(detector, position, x, y, fraction);

	  if (x[0] != INVALID_PIXEL) {
	    struct Event event;
	  
	    // Store the photon charge and the new arrival time:
	    event.pha = get_pha(energy, detector);  // TODO: RMF
	    event.time = time;
	    event.xi = x[0];
	    event.yi = y[0];
	    event.grade = 0;
	    event.frame = detector->frame;

	    // Add the event to the FITS event list.
	    if ((event.pha>=detector->pha_threshold)&&
		(energy>=detector->energy_threshold)){ // Check lower PHA threshold
	      // There is an event in this pixel, so insert it into eventlist:
	      add_eventlist_row(&eventlist_file, event, &status);
	    }
	  } // END x[0] != INVALID_PIXEL

	} // END detector->type == TES

      } // END 'time' within specified time interval

    }  // END of scanning the impact list.

  } while(0);  // END of the error handling loop.

  // --- END of Detection process ---


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Release HEADAS random number generator
  HDmtFree();

  if (impactlist_fptr) fits_close_file(impactlist_fptr, &status);
  if (eventlist_file.fptr) fits_close_file(eventlist_file.fptr, &status);

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

  // Release memory of detector Redistribution Matrix
  free_rmf(&detector->rmf);
  
  // Release memory of detector EBOUNDS
  free_ebounds(&detector->ebounds);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully\n\n");

  return(status);
}




////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int photon_detection_getpar(
		       char impactlist_filename[],
		       char rmf_filename[],       
		       char eventlist_filename[],
		       double *t0,        // start time
		       double *timespan,  // time span
		       Detector *detector
		       )
{
  char msg[MAXMSG];             // error output buffer
  int status=EXIT_SUCCESS;      // error status

  // Get the name of the impact list file (FITS file)
  if ((status = PILGetFname("impactlist_filename", impactlist_filename))) {
    sprintf(msg, "Error reading the name of the impact list file!\n");
    HD_ERROR_THROW(msg,status);
    return(status);
  }

  // Get the detector type (FRAMESTORE, DEPFET, TES microcalorimeter, ...)
  int type;
  if ((status = PILGetInt("detectortype", &type))) {
    sprintf(msg, "Error reading the detector type!\n");
    HD_ERROR_THROW(msg,status);
    return (status);

  } else {
    switch (type) {
      // According to the different detector types, the program
      // has to read different parameters.
    case 1: 
      detector->type = FRAMESTORE; 

      // Get the integration time for the FRAMESTORE CCD.
      if ((status = PILGetReal("integration_time", &detector->integration_time))) {
	sprintf(msg, "Error reading the integration time!\n");
	HD_ERROR_THROW(msg,status);
      }

      break;

    case 2: 
      detector->type = DEPFET;     

      // Get the dead time for the DEPFET APS (readout time per line).
      if ((status = PILGetReal("dead_time", &detector->dead_time))) {
	sprintf(msg, "Error reading the dead time!\n");
	HD_ERROR_THROW(msg,status);
      } 
      // Get the clear time for the DEPFET APS (time required to clear one line).
      else if ((status = PILGetReal("clear_time", &detector->clear_time))) {
	sprintf(msg, "Error reading the clear time!\n");
	HD_ERROR_THROW(msg,status);
      }
      
      break;

    case 3: 
      detector->type = TES;        
      break;

    case 4:
      detector->type = HTRS;

      // Get the dead time for the HTRS (readout time of the charge 
      // cloud in a pixel).
      if ((status = PILGetReal("dead_time", &detector->dead_time))) {
	sprintf(msg, "Error reading the dead time!\n");
	HD_ERROR_THROW(msg,status);
      } 

      break;

    default:     
      detector->type = 0;
      sprintf(msg, "Error: incorrect detector type!\n");
      HD_ERROR_THROW(msg,status);
    }
  }
  if (status) return(status);
  // END of handling different detector types.


  // detector width [pixel]
  if (detector->type == HTRS) {
    detector->width = 7;  // HTRS can only be used for fixed pixel size!
  } else {
    if ((status = PILGetInt("det_width", &detector->width))) {
      sprintf(msg, "Error reading the width of the detector!\n");
      HD_ERROR_THROW(msg,status);
      return(status);
    }
  }

  // [m]
  if ((status = PILGetReal("det_pixelwidth", &detector->pixelwidth))) {
    sprintf(msg, "Error reading the width of the detector pixels!\n");
    HD_ERROR_THROW(msg,status);
  }

  // [m]
  else if ((status = PILGetReal("ccsigma", &detector->ccsigma))) {
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
    detector->pha_threshold = (long)pha_threshold;
  }
  if (detector->pha_threshold==0) {
    if ((status = PILGetReal4("energy_threshold", &detector->energy_threshold))) {
      sprintf(msg, "Error: could not determine detector energy threshold!\n");
      HD_ERROR_THROW(msg,status);
      return(status);
    }
  } else {
    detector->energy_threshold=0.;
  }

  // Get the name of the detector redistribution file (FITS file)
  if ((status = PILGetFname("rmf_filename", rmf_filename))) {
    sprintf(msg, "Error reading the name of the detector" 
	    "redistribution matrix file (RMF)!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the name of the output event list (FITS file)
  else if ((status = PILGetFname("eventlist_filename", eventlist_filename))) {
    sprintf(msg, "Error reading the name of the event list file!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the start time of the simulation
  else if ((status = PILGetReal("t0", t0))) {
    sprintf(msg, "Error reading the 't0' parameter!\n");
    HD_ERROR_THROW(msg,status);
  }

  // Get the timespan for the simulation
  else if ((status = PILGetReal("timespan", timespan))) {
    sprintf(msg, "Error reading the 'timespan' parameter!\n");
    HD_ERROR_THROW(msg,status);
  }


  return(status);
}


