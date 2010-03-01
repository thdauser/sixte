#include "framestoredetector.h"


int initFramestoreDetector(FramestoreDetector* fd, 
			   struct FramestoreDetectorParameters* parameters)
{
  int status = EXIT_SUCCESS;

  // Call the initialization routines of the underlying data structures.
  // Init the generic detector properties like the RMF.
  status = initGenericDetector(&fd->generic, &parameters->generic);
  if (EXIT_SUCCESS!=status) return(status);
#ifdef EXPONENTIAL_SPLITS
  headas_chat(5, "exponential split model\n");
#else
  headas_chat(5, "Gaussian split model\n");
#endif
  // Get the memory for the pixel array.
  status = initSquarePixels(&fd->pixels, &parameters->pixels);
  if (EXIT_SUCCESS!=status) return(status);
  
  // Set up the framestore configuration.
  fd->integration_time = parameters->integration_time;

  // Set the first readout time such that the first readout is performed 
  // immediately at the beginning of the simulation.
  fd->readout_time = parameters->t0;
  fd->frame = 0;

  // Create and open new event list file.
  status = openNeweROSITAEventFile(&fd->eventlist, parameters->eventlist_filename,
				       parameters->eventlist_template);
  if (EXIT_SUCCESS!=status) return(status);

  return(status);
}



int cleanupFramestoreDetector(FramestoreDetector* fd) 
{
  int status=EXIT_SUCCESS;

  // Call the cleanup routines of the underlying data structures.
  cleanupSquarePixels(&fd->pixels);
  cleanupGenericDetector(&fd->generic);
  status+=closeeROSITAEventFile(&fd->eventlist);

  return(status);
}



int checkReadoutFramestoreDetector(FramestoreDetector* fd, double time)
{
  int status = EXIT_SUCCESS;

  // Check, if the detector integration time was exceeded. 
  // In that case, read out the detector.
  if (time > fd->readout_time) {
    // It's time to read out the detector!
      
    // Readout the detector and create eventlist entries for the actual time:
    status = readoutFramestoreDetector(fd);
    if (EXIT_SUCCESS!=status) return(status);

    // Clear the detector array.
    clearSquarePixels(&fd->pixels);

    // Update the detector frame time to the next frame until the current
    // time is within the detector->readout interval.
    // This CAN ONLY BE DONE for FRAMESTORE detectors!
    // For detectors with individual readout lines a more complicated method is required.
    while (time > fd->readout_time) {
      fd->readout_time += fd->integration_time;
      fd->frame += 2;
    }

    // Print the time of the current events (program status
    // information for the user).
    headas_printf("\rtime: %.3lf s ", fd->readout_time);
    fflush(NULL);
  }

  return(status);
}



inline int readoutFramestoreDetector(FramestoreDetector* fd) 
{
  int x, y;          // Counters for the loop over the pixel array.

  eROSITAEvent list[100]; // List of events belonging to the same pattern.
  int nlist=0, count;     // Number of entries in the list.

  int status = EXIT_SUCCESS;


  // Find the events in the pixel array.
  for (x=0; x<fd->pixels.xwidth; x++) {
    for (y=0; y<fd->pixels.ywidth; y++) {
      
      // Check if the pixel contains any charge. If there is no charge in it at
      // all, there is no use to determine the PHA channel.
      if (fd->pixels.array[x][y].charge > 1.e-6) {

	// Determine the detector channel that corresponds to the charge stored
	// in the detector pixel.
	list[0].pha = getChannel(fd->pixels.array[x][y].charge, fd->generic.rmf);
	
	// The PHA channel should only be less than zero, when the photon 
	// is lost, i.e., not detected at all. As the RSP is usually normalized,
	// i.e., it only contains the RMF, this should never be the case.
	assert(list[0].pha >= 0);
	// Maybe: if (event.pha < 0) continue;
	
	// Check lower threshold (PHA and energy):
	if ((list[0].pha>=fd->generic.pha_threshold) && 
	    (fd->pixels.array[x][y].charge>=fd->generic.energy_threshold)) { 
	  
	  // The event is the first in the list (containing all events
	  // that belong to the same pattern).
	  nlist = 1;
	  list[0].energy = fd->pixels.array[x][y].charge * 1.e3; // [eV]
	  list[0].xi = x;
	  list[0].yi = y;
	  // Delete the charge of this event.
	  fd->pixels.array[x][y].charge = 0.;

	  // Call marker routine to check for surrounding pixels.
	  fdMarkEvents(list, &nlist, fd, x, y);
	  
	  // Perform pattern type identification.
	  for (count=0; count<nlist; count++) {
	    // TODO
	    list[count].pat_typ = nlist;
	  }

	  // Store the list in the event file.
	  for (count=0; count<nlist; count++) {
	    // Set missing properties.
	    list[count].time  = fd->readout_time;
	    list[count].frame = fd->frame;
	    
	    status=addeROSITAEvent2File(&fd->eventlist, &(list[count]));
	    if (EXIT_SUCCESS!=status) return(status);
	  }
	  
	} // End of check if event is above specified threshold.
      } // END of check if pixel contains any charge.
    } // END of loop over x
  } // END of loop over y
  // END of find split patterns.

  return(status);
}



int addImpact2FramestoreDetector(FramestoreDetector* fd, Impact* impact)
{
  int status=EXIT_SUCCESS;

  // Before adding the new impact to the detector check whether
  // a readout has to be performed in advance.
  status=checkReadoutFramestoreDetector(fd, impact->time);

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
    
    // Determine the affected detector pixels (including split partners).
#ifdef EXPONENTIAL_SPLITS
    int npixels = getSquarePixelsExponentialSplits(&fd->pixels, &(fd->generic.ecc), 
						   impact->position, x, y, fraction);
#else
    int npixels = getSquarePixelsGaussianSplits(&fd->pixels, &(fd->generic.gcc), 
						impact->position, x, y, fraction);
#endif
    
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

  return(status);
}



/*
// Consistency check for size of charge cloud:
if (detector->ccsize > detector->pixelwidth) {
status=EXIT_FAILURE;
HD_ERROR_THROW("Error: charge cloud size greater than pixel width!\n", status);
break;
}
*/



void fdMarkEvents(eROSITAEvent* list, int* nlist, 
		  FramestoreDetector* fd, 
		  int x, int y)
{
  // Split threshold. TODO
  const double split_threshold = 1.e-2;

  // Possible neighbors (no diagonal neighbors).
  const int neighbors_x[4] = {0, 0, 1, -1};
  const int neighbors_y[4] = {1, -1, 0, 0};

  int neighbor;

  for (neighbor=0; neighbor<4; neighbor++) {
    assert(*nlist+1 < 100);

    // Check if the neighbor is within the detector dimensions.
    if ((x+neighbors_x[neighbor]<0) || 
	(x+neighbors_x[neighbor]>=fd->pixels.xwidth) ||
	(x+neighbors_x[neighbor]<0) || 
	(x+neighbors_x[neighbor]>=fd->pixels.xwidth))
      continue;

    if (fd->pixels.array[x+neighbors_x[neighbor]][y+neighbors_y[neighbor]].charge > 
	split_threshold) {

      (*nlist)++;
      printf("nlist: %d (%d, %d)\n", *nlist, 
	     x+neighbors_x[neighbor],
	     y+neighbors_y[neighbor]);

      list[*nlist-1].pha = 
	getChannel(fd->pixels.array[x+neighbors_x[neighbor]][y+neighbors_y[neighbor]].charge, 
		   fd->generic.rmf);
      list[*nlist-1].energy = 
	fd->pixels.array[x+neighbors_x[neighbor]][y+neighbors_y[neighbor]].charge 
	* 1.e3; // [eV]
      list[*nlist-1].xi = x+neighbors_x[neighbor];
      list[*nlist-1].yi = y+neighbors_y[neighbor];

      // Clear the charge from the array.
      fd->pixels.array[x+neighbors_x[neighbor]][y+neighbors_y[neighbor]].charge = 0.;

      // Call marker routine to check for surrounding pixels.
      fdMarkEvents(list, nlist, fd, 
		   x+neighbors_x[neighbor], 
		   y+neighbors_y[neighbor]);

    } // END of check if event is above the split threshold.
  } // END of loop over all neighbors.

}
