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
  headas_chat(5, "exponential model (by Konrad Dennerl) for split events\n");
#else
  headas_chat(5, "Gaussian charge cloud model for split events\n");
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
  int x, y;               // Counters for the loop over the pixel array.
  long pha;               // Buffer for PHA values.

  eROSITAEvent list[100]; // List of events belonging to the same pattern.
  int nlist=0, count;     // Number of entries in the list.

  int maxidx, minidx;     // Indices of the maximum and minimum event in the 
                          // event list.

  int status = EXIT_SUCCESS;


  // Find the events in the pixel array.
  for (x=0; x<fd->pixels.xwidth; x++) {
    for (y=0; y<fd->pixels.ywidth; y++) {
      
      // Check if the pixel contains any charge. If there is no charge in it at
      // all, there is no use to determine the PHA channel.
      if (fd->pixels.array[x][y].charge > 1.e-6) {

	// Determine the detector channel that corresponds to the charge stored
	// in the detector pixel.
	pha = getChannel(fd->pixels.array[x][y].charge, fd->generic.rmf);
	
	// The PHA channel should only be less than zero, when the photon 
	// is lost, i.e., not detected at all. As the RSP is usually normalized,
	// i.e., it only contains the RMF, this should never be the case.
	assert(pha >= 0);
	// Maybe: if (pha < 0) continue;
	
	// Check lower threshold (PHA and energy):
	if ((pha>=fd->generic.pha_threshold) && 
	    (fd->pixels.array[x][y].charge>=fd->generic.energy_threshold)) { 

	  // Clear the event list.
	  nlist=0;
	  maxidx=0;
	  minidx=0;
	  
	  // Call marker routine to check for surrounding pixels.
	  fdMarkEvents(list, &nlist, &maxidx, &minidx, fd, x, y);

	  // Perform pattern type identification.
	  fdPatternIdentification(list, nlist, maxidx, minidx);

	  // Store the list in the event FITS file.
	  for (count=0; count<nlist; count++) {
	    // Set missing properties.
	    list[count].time  = fd->readout_time;
	    list[count].frame = fd->frame;
	    
	    status=addeROSITAEvent2File(&fd->eventlist, &(list[count]));
	    if (EXIT_SUCCESS!=status) return(status);

	  } // End of storing the event list in the FITS file.
	  
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
		  int* maxidx, int* minidx,
		  FramestoreDetector* fd, 
		  int x, int y)
{
  // Split threshold. TODO
  const double split_threshold = 1.e-2;

  // Possible neighbors (no diagonal neighbors).
  const int neighbors_x[4] = {0, 0, 1, -1};
  const int neighbors_y[4] = {1, -1, 0, 0};

  int neighbor, xi, yi;

  // Create a new event in the list.
  assert(*nlist+1 < 100);  
  list[*nlist].pha = getChannel(fd->pixels.array[x][y].charge, fd->generic.rmf);
  list[*nlist].energy = fd->pixels.array[x][y].charge * 1.e3; // [eV]
  list[*nlist].xi  = x;
  list[*nlist].yi  = y;
  (*nlist)++;
  // Delete the charge of this event from the pixel array.
  fd->pixels.array[x][y].charge = 0.;

  // Check if the new event has the maximum or minium energy in the event list.
  if (list[*nlist-1].energy <= list[*minidx].energy) {
    *minidx = *nlist-1;
  } else if (list[*nlist-1].energy > list[*maxidx].energy) {
    *maxidx = *nlist-1;
  }

#ifdef FD_DETECT_PATTERNS
  // Loop over the directly neighboring pixels.
  for (neighbor=0; neighbor<4; neighbor++) {

    // Coordinates of the neighboring pixel:
    xi = x+neighbors_x[neighbor];
    yi = y+neighbors_y[neighbor];

    // Check if the neighbor is within the detector dimensions.
    if ((xi<0) || (xi>=fd->pixels.xwidth) ||
	(yi<0) || (yi>=fd->pixels.ywidth))
      continue;

    if (fd->pixels.array[xi][yi].charge > split_threshold) {

      // Call marker routine to check for surrounding pixels.
      fdMarkEvents(list, nlist, maxidx, minidx, fd, xi, yi);

    } // END of check if event is above the split threshold.

  } // END of loop over all neighbors.
#endif

}



void fdPatternIdentification(eROSITAEvent* list, int nlist, 
			     int maxidx, int minidx)
{
 
  // Single events.
  if (1==nlist) {
    list[0].pat_inf = 0;
  }

  // Double events.
  else if (2==nlist) {
    // 0i0
    // 0m0
    // 000
    if ((list[maxidx].xi==list[minidx].xi) && 
	(list[maxidx].yi==list[minidx].yi-1)) {
      list[maxidx].pat_inf = 15;
      list[minidx].pat_inf = 12;
    }
    // 000
    // 0mi
    // 000
    else if ((list[maxidx].xi==list[minidx].xi-1) && 
	     (list[maxidx].yi==list[minidx].yi)) {
      list[maxidx].pat_inf = 15;
      list[minidx].pat_inf = 16;
    }
    // 000
    // 0m0
    // 0i0
    else if ((list[maxidx].xi==list[minidx].xi) && 
	     (list[maxidx].yi==list[minidx].yi+1)) {
      list[maxidx].pat_inf = 15;
      list[minidx].pat_inf = 18;
    }
    // 000
    // im0
    // 000
    else if ((list[maxidx].xi==list[minidx].xi+1) && 
	     (list[maxidx].yi==list[minidx].yi)) {
      list[maxidx].pat_inf = 15;
      list[minidx].pat_inf = 14;
    }
    // Invalid pattern. (Actually this code should never be executed.)
    else {
      printf("Invalid double event!\n");
      list[maxidx].pat_inf = 0;
      list[minidx].pat_inf = 0;
    }
  } // END of double events.
  
}
