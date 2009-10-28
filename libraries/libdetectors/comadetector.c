#include "comadetector.h"


CoMaDetector* initCoMaDetector(struct CoMaDetectorParameters* parameters,
			       int* status)
{
  // Allocate memory for the CoMaDetector object.
  CoMaDetector* det=(CoMaDetector*)malloc(sizeof(CoMaDetector));
  if (NULL==det) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for CoMaDetector "
		   "object!\n", *status);
    return(det);
  }
  // Set initial values.
  det->pixels    = NULL;
  det->eventfile = NULL;

  // Call the initialization routines of the underlying data structures.
  det->pixels = getSquarePixels(&parameters->pixels, status);
  if (EXIT_SUCCESS!=*status) return(det);

  // Create and open new event list file.
  det->eventfile = openNewCoMaEventFile(parameters->eventfile_filename,
					parameters->eventfile_template,
					status);
  if (EXIT_SUCCESS!=*status) return(det);

  return(det);
}



void freeCoMaDetector(CoMaDetector* det) 
{
  if (NULL!=det) {
    // Call the clean-up routines of the underlying data structures.
    freeSquarePixels(det->pixels);
    closeCoMaEventFile(det->eventfile);
    // Free the allocated memory.
    free(det);
  }
}



int addImpact2CoMaDetector(CoMaDetector* det, Impact* impact)
{
  int status=EXIT_SUCCESS;

  /*
  // Before adding the new impact to the detector check whether
  // a readout has to be performed in advance.
  status=checkReadoutCoMaDetector(fd, impact->time);

  // Determine a detector channel (PHA channel) according to the RMF.
  // The channel is obtained from the RMF using the corresponding
  // HEAdas routine which is based on drawing a random number.
  long channel;
  ReturnChannel(fd->generic.rmf, impact->energy, 1, &channel);

  // Check if the photon is really measured. If the
  // PHA channel returned by the HEAdas RMF function is '-1', 
  // the photon is not detected.
  // This can happen as the RMF is actually an RSP including, e.g., 
  // the detector quantum efficiency and filter transmission.
  if (-1==channel) {
    return(status); // Break the function (and continue with the next photon).
  }
  assert(channel>=0);

  // Get the corresponding created charge.
  // NOTE: In this simulation the charge is represented by the nominal
  // photon energy which corresponds to the PHA channel according to the
  // EBOUNDS table.
  float charge = getEnergy(channel, fd->generic.rmf);
  
  if (charge > 0.) {
    int x[4], y[4];
    double fraction[4];
    
    // Determine the affected detector pixels.
    int npixels = getSquarePixelsSplits(&fd->pixels, &fd->generic, impact->position, 
					x, y, fraction);
    
    // Add the charge created by the photon to the affected detector pixels.
    int count;
    for (count=0; count<npixels; count++) {
      if (x[count] != INVALID_PIXEL) {
	fd->pixels.array[x[count]][y[count]].charge += 
	  charge * fraction[count];
	  // |      |-> charge fraction due to split events
	  // |-> charge created by incident photon
      }
    }
  } // END if(charge>0.)
  */

  return(status);
}


