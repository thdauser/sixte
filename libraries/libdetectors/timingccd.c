#include "timingccd.h"


int initTimingccd(TimingCCD* tc, struct TimingCCDParameters* parameters)
{
  int status = EXIT_SUCCESS;

  // Set up the initial CCD configuration.
  tc->step_time = parameters->step_time;
  // Set the first readout time such that the first readout is performed 
  // immediately at the beginning of the simulation.
  tc->time = parameters->t0;
  
  // Call the initialization routines of the underlying data structures.
  // Init the generic detector properties like the RMF.
  status = initGenericDetector(&tc->generic, &parameters->generic);
  if (EXIT_SUCCESS!=status) return(status);
  // Get the memory for the pixel array.
  status = initSquarePixels(&tc->pixels, &parameters->pixels);
  if (EXIT_SUCCESS!=status) return(status);
  
  return(status);
}



int cleanupTimingCCD(TimingCCD* tc) 
{
  int status=EXIT_SUCCESS;

  // Call the cleanup routines of the underlying data structures.
  cleanupSquarePixels(&tc->pixels);
  cleanupGenericDetector(&tc->generic);

  return(status);
}


/*
int checkReadoutTimingccd(Timingccd* fd, double time)
{
  int status = EXIT_SUCCESS;

  // Check, if the detector integration time was exceeded. 
  // In that case, read out the detector.
  if (time > fd->readout_time) {
    // It's time to read out the detector!
      
    // Readout the detector and create eventlist entries for the actual time:
    status = readoutTimingccd(fd);
    if (EXIT_SUCCESS!=status) return(status);

    // Clear the detector array.
    clearSquarePixels(&fd->pixels);

    // Update the detector frame time to the next frame until the current
    // time is within the detector->readout interval.
    // This CAN ONLY BE DONE for FRAMESTORE detectors!
    // For detectors with individual readout lines a more complicated method is required.
    while (time > fd->readout_time) {
      fd->readout_time += fd->integration_time;
      fd->frame += 1;
    }

    // Print the time of the current events (program status
    // information for the user).
    headas_printf("\rtime: %.3lf s ", fd->readout_time);
    fflush(NULL);
  }

  return(status);
}


float fdMaximumCharge(float* charges, int nlist)
{
  int count;
  float maximum=0.;
  for(count=0; count<nlist; count++) {
    if (charges[count]>maximum) {
      maximum = charges[count];
    }
  }
  return(maximum);
}



inline int readoutTimingccd(Timingccd* fd) 
{
  int x, y;   // Counters for the loop over the pixel array.
  int x2, y2;
  long pha;   // Buffer for PHA values.

  // List of events belonging to the same pattern.
  eROSITAEvent list[MAX_N_SPLIT_LIST]; 
  int nlist=0, count; // Number of entries in the list.
  int maxidx, minidx; // Indices of the maximum and minimum event in the 
                      // event list.

  int status = EXIT_SUCCESS;


  // Find the events in the pixel array.
  for (x=0; x<fd->pixels.xwidth; x++) {
    for (y=0; y<fd->pixels.ywidth; y++) {
      
      // Check if the pixel contains any charge. If there is no charge 
      // in it at all, there is no use to determine the PHA channel.
      if (fd->pixels.array[x][y].charge > 1.e-6) {

	// Determine the detector channel that corresponds to the 
	// charge stored in the detector pixel.
	pha = getChannel(fd->pixels.array[x][y].charge, fd->generic.rmf);
	
	// The PHA channel should only be less than zero, when the photon 
	// is lost, i.e., not detected at all. As the RSP is usually 
	// normalized,
	// i.e., it only contains the RMF, this should never be the case.
	assert(pha >= 0);
	// Maybe: if (pha < 0) continue;
	
	// Check lower threshold (PHA and energy):
	if ((pha>=fd->generic.pha_threshold) && 
	    (fd->pixels.array[x][y].charge>=fd->generic.energy_threshold)) { 

	  // Add the central pixel to the list.
	  nlist = 1;
	  maxidx= 0;
	  minidx= 0;
	  list[0].xi = x;
	  list[0].yi = y;
	  list[0].energy = fd->pixels.array[x][y].charge * 1.e3; // [eV]
	  list[0].pha    = pha;

	  // Delete the pixel charge after it has been read out.
	  fd->pixels.array[x][y].charge = 0.;

	  // Check the surrounding pixels (according to the method 
	  // proposed by Konrad Dennerl for the on-board processor).
	  float neighbor_charges[4] = { 0., 0., 0., 0. };
	  float combined_neighbor_charges[4] = { 0., 0., 0., 0. };
	  // Right:
	  if (x+1<fd->pixels.xwidth) {
	    neighbor_charges[0] = fd->pixels.array[x+1][y].charge;
	  }
	  // Left:
	  if (x-1>=0) {
	    neighbor_charges[1] = fd->pixels.array[x-1][y].charge;
	  }
	  // Top:
	  if (y+1<fd->pixels.ywidth) {
	    neighbor_charges[2] = fd->pixels.array[x][y+1].charge;
	  }
	  // Bottom:
	  if (y-1>=0) {
	    neighbor_charges[3] = fd->pixels.array[x][y-1].charge;
	  }
	  // Calculate the quantities c_rt, c_tl, c_lb, and c_br:
	  // c_rt:
	  combined_neighbor_charges[0] = 
	    neighbor_charges[0] + neighbor_charges[2];
	  // c_tl:
	  combined_neighbor_charges[1] = 
	    neighbor_charges[2] + neighbor_charges[1];
	  // c_lb:
	  combined_neighbor_charges[2] = 
	    neighbor_charges[1] + neighbor_charges[3];
	  // c_br:
	  combined_neighbor_charges[3] = 
	    neighbor_charges[3] + neighbor_charges[0];

	  // Determine the split threshold as the sum of the charge of the 
	  // central pixel and the maximum surrounding charge times the 
	  // selected threshold fraction.
	  float split_threshold = 
	    (fdMaximumCharge(combined_neighbor_charges, 4) + list[0].energy*1.e-3) *
	    //                                               [eV] -> [keV]
	    fd->split_threshold;

	  // Check the 3x3 environment for pixels with a charge above the 
	  // split threshold.
	  for (x2=MAX(x-1,0); x2<MIN(x+2,fd->pixels.xwidth); x2++) {
	    for (y2=MAX(y-1,0); y2<MIN(y+2, fd->pixels.ywidth); y2++) {

	      // Do not regard the central pixel. This has already been 
	      // added to the list.
	      if ((x==x2)&&(y==y2)) continue;

	      if (fd->pixels.array[x2][y2].charge > split_threshold) {
		list[nlist].xi = x2;
		list[nlist].yi = y2;
		list[nlist].energy = 
		  fd->pixels.array[x2][y2].charge * 1.e3; // [eV]
		list[nlist].pha = 
		  getChannel(fd->pixels.array[x2][y2].charge, fd->generic.rmf);

		// Delete the pixel charge after it has been read out.
		fd->pixels.array[x2][y2].charge = 0.;

		// Check if the new event has the maximum or minium 
		// energy in the split list.
		if (list[nlist].energy <= list[minidx].energy) {
		  minidx = nlist;
		} else if (list[nlist].energy > list[maxidx].energy) {
		  maxidx = nlist;
		}

		nlist++;
	      } // TODO:
	      // else { 
	      // delete pixel charge (otherwise it might happen, that split
	      // partners lie below the split threshold, but are above the 
	      // event threshold and are therefore read out in the next step
	      // resulting in an invalid pattern (double or triple).
	      // }
	    }
	  }

	  // Perform the pattern type identification.
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
*/


int addImpact2Timingccd(TimingCCD* tc, Impact* impact)
{
  int status=EXIT_SUCCESS;

  // Before adding the new impact to the detector check whether
  // a readout has to be performed in advance.
  // TODO
//  status=checkReadoutTimingCCD(tc, impact->time);

  // Determine a detector channel (PHA channel) according to the RMF.
  // The channel is obtained from the RMF using the corresponding
  // HEAdas routine which is based on drawing a random number.
  long channel;
  ReturnChannel(tc->generic.rmf, impact->energy, 1, &channel);

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
  float charge = getEnergy(channel, tc->generic.rmf);
  
  if (charge > 0.) {
    int x[4], y[4];
    double fraction[4];
    
    // Determine the affected detector pixels (including split partners).
#ifdef TIMING_CCD_EXPONENTIAL_SPLITS
    int npixels = getSquarePixelsExponentialSplits(&tc->pixels, &(tc->generic.ecc), 
						   impact->position, x, y, fraction);
#else
    int npixels = getSquarePixelsGaussianSplits(&tc->pixels, &(tc->generic.gcc), 
						impact->position, x, y, fraction);
#endif
    
    // Add the charge created by the photon to the affected detector pixels.
    int count;
    for (count=0; count<npixels; count++) {
      if (x[count] != INVALID_PIXEL) {
	tc->pixels.array[x[count]][y[count]].charge += 
	  charge * fraction[count];
	  // |      |-> charge fraction due to split events
	  // |-> charge created by incident photon
      }
    }
  } // END if(charge>0.)

  return(status);
}



