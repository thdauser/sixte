#include "framestore.h"


int initFramestoreDetector(FramestoreDetector* fd, 
			   struct FramestoreDetectorParameters* parameters)
{
  int status = EXIT_SUCCESS;

  // Call the initialization routines of the underlying data structures.
  initGenericDetector(&fd->generic, &parameters->generic);
  initSquarePixels(&fd->pixels, &parameters->pixels);

  // Set up the framestore configuration.
  fd->integration_time = parameters->integration_time;

  // Set the first readout time such that the first readout is performed 
  // immediately at the beginning of the simulation.
  fd->readout_time = parameters->t0;
  fd->frame = 0;

  return(status);
}



void readoutFramestoreDetector(FramestoreDetector* fd) 
{
  return;
}

void addImpact2FramestoreDetector(FramestoreDetector* fd)
{
  return;
}


/*
// Consistency check for size of charge cloud:
if (detector->ccsize > detector->pixelwidth) {
status=EXIT_FAILURE;
HD_ERROR_THROW("Error: charge cloud size greater than pixel width!\n", status);
break;
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

*/

