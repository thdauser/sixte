#include "framestore.h"


//////////////////////////////////////////////////
int init_FramestoreDetector(Detector* detector, struct FramestoreParameters parameters) {
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

    // Set the framestore-specific properties.
    properties->frame = 0;
    properties->integration_time = parameters.integration_time;

    // Set the charge cloud size:
    detector->ccsigma = parameters.ccsigma;
    detector->ccsize  = 3.*parameters.ccsigma;

    // Consistency check for size of charge cloud:
    if (detector->ccsize > detector->pixelwidth) {
      status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: charge cloud size greater than pixel width!\n", status);
      break;
    }

    // Set the first readout time such that the first readout is performed 
    // immediately at the beginning of the simulation.
    detector->readout_time = parameters.t0;

    // Set the readout routine:
    detector->readout = readout_FramestoreDetector;

    // Get the memory for the detector pixels
    if (get_DetectorPixels(detector, &status)) break;

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
