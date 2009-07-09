#include "framestore.h"


//////////////////////////////////////////////////
int init_FramestoreDetector(Detector* detector, struct DetectorParameters detpar,
			    struct FramestoreParameters framepar) 
{
  struct FramestoreProperties* properties = NULL;
  int status = EXIT_SUCCESS;
  
  do { // Error handling loop.

    // Set the detector tpye.
    detector->type = FRAMESTORE;

    // Allocate memory for the framestore-specific elements in the detector data
    // structure.
    properties = malloc(sizeof(struct FramestoreProperties));
    if (NULL==properties) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for detector specific elements failed !\n", 
		     status);
      break;
    }

    // Set the GENERAL detector properties.
    // Detector size:
    detector->width  = detpar.width;
    detector->offset = detpar.width/2;
    detector->pixelwidth = detpar.pixelwidth;

    // Set the charge cloud size:
    detector->ccsigma =    detpar.ccsigma;
    detector->ccsize  = 3.*detpar.ccsigma;

    // Consistency check for size of charge cloud:
    if (detector->ccsize > detector->pixelwidth) {
      status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: charge cloud size greater than pixel width!\n", status);
      break;
    }
    
    detector->frame = 0;

    // Thresholds:
    detector->pha_threshold = detpar.pha_threshold;
    detector->energy_threshold = detpar.energy_threshold;

    // Get the memory for the detector pixels
    if (get_DetectorPixels(detector, &status)) break;

    // Read the detector RMF and EBOUNDS from the specified file and 
    // assign them to the Detector data structure.
    if ((status=detector_assign_rsp(detector, detpar.rmf_filename)) 
	!= EXIT_SUCCESS) break;


    // Set the FRAMESTORE-SPECIFIC properties.
    //    properties->frame = 0;
    properties->integration_time = framepar.integration_time;

    // Set the first readout time such that the first readout is performed 
    // immediately at the beginning of the simulation.
    detector->readout_time = detpar.t0;


    // Set the readout routine:
    detector->readout = readout_FramestoreDetector;

    // Set the photon detection routine:
    detector->add_impact = add_Impact2FramestoreDetector;

  } while(0); // End of Error handling loop.

  detector->specific = properties; // (void*)
  return(status);
}




//////////////////////////////////////////////////////////////////////
void readout_FramestoreDetector(
				void* det,
				double time, 
				struct Eventlist_File* eventlist_file,
				int *status
				) 
{
  Detector* detector = (Detector*) det;

  // Check, if the detector integration time was exceeded. 
  // In that case, read out the detector.
  if (time > detector->readout_time) {

    // Readout the detector and create eventlist entries for the actual time:
    readout(detector, eventlist_file, status);

    // Clear the detector array.
    clear_detector(detector);

    // Update the detector frame time to the next frame until the current
    // time is within the detector->readout interval.
    // This CAN ONLY BE DONE for FRAMESTORE detectors!!
    // For detectors with individual readout lines a more complicated method is required.
    struct FramestoreProperties* specific = (struct FramestoreProperties*)detector->specific;
    while (time > detector->readout_time) {
      detector->readout_time += specific->integration_time;
      detector->frame+=2;
    }

    // Print the time of the current events in order (program status
    // information for the user).
    headas_chat(0, "\rtime: %.3lf s ", detector->readout_time);
    fflush(NULL);
  }

}




/////////////////////////////////////////////////////
void add_Impact2FramestoreDetector (void* det, struct Impact* impact) {
  Detector* detector = (Detector*) det;

  // Determine a detector channel (PHA channel) according to RMF.
  // The channel is obtained from the RMF using the corresponding
  // HEAdas routine which is based on drawing a random number.
  long channel;
  ReturnChannel(detector->rmf, impact->energy, 1, &channel);

  // Check if the photon is really measured. If the
  // PHA channel returned by the HEAdas RMF function is '-1', 
  // the photon is not detected.
  if (channel==-1) {
    return; // Break the function (and continue with the next photon).
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
    int npixels = get_pixel_square(detector, impact->position, x, y, fraction);
    
    // Add the charge created by the photon to the affected detector pixels.
    int count;
    for (count=0; count<npixels; count++) {
      if (x[count] != INVALID_PIXEL) {
	detector->pixel[x[count]][y[count]].charge += 
	  charge * fraction[count] * 
	  // |      |-> charge fraction due to split events
	  // |-> charge created by incident photon
	  detector_active(x[count], y[count], detector, impact->time);
	// |-> "1" if pixel can measure charge, "0" else
      }
    }
  } // END if(charge>0.)

}



