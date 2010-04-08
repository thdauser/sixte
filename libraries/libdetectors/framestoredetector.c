#include "framestoredetector.h"


int initFramestoreDetector(FramestoreDetector* fd, 
			   struct FramestoreDetectorParameters* parameters)
{
  int status = EXIT_SUCCESS;

  // Set up the framestore configuration.
  fd->integration_time = parameters->integration_time;
  fd->split_threshold  = parameters->split_threshold;
  // Set the first readout time such that the first readout is performed 
  // immediately at the beginning of the simulation.
  fd->readout_time = parameters->t0;
  fd->frame = 0;
  
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
  
  // Create and open new event list file.
  status = openNeweROSITAEventFile(&fd->eventlist, 
				   parameters->eventlist_filename,
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


/* Deprecated 
float fdGetPixelCharge(FramestoreDetector* fd, int x, int y) 
{
  // Check if the pixel coordinates lie within the valid array.
  if ((x<0) || (x >= fd->pixels.xwidth) || 
      (y<0) || (y >= fd->pixels.ywidth)) {
    return (0.);
  } else {
    return (fd->pixels.array[x][y].charge);
  }
}
*/


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



inline int readoutFramestoreDetector(FramestoreDetector* fd) 
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

	  /* Deprecated routine to detect split patterns. The routine
	     checks not only the surrounding 3x3 environment of pixels
	     but is recursively called in order to identify larger invalid
	     patterns. 
	  // Clear the event list.
	  nlist=0;
	  maxidx=0;
	  minidx=0;
	  
	  // Call marker routine to check the surrounding pixels for 
	  // events above the split threshold.
	  fdMarkEvents(list, &nlist, &maxidx, &minidx, fd, x, y);

	  // TODO Determine split threshold from charge of MAXIMUM split partner.
	  // Maybe insert a new routine here!
	  printf("MarkEvents: nlist = %d\n", nlist);

	  // Perform pattern type identification.
	  fdPatternIdentification(list, nlist, maxidx, minidx);
	  */
	  
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
	  /*
	    fdGetPixelCharge(fd, x+1, y) + fdGetPixelCharge(fd, x, y+1),
	    fdGetPixelCharge(fd, x-1, y) + fdGetPixelCharge(fd, x, y+1),
	    fdGetPixelCharge(fd, x-1, y) + fdGetPixelCharge(fd, x, y-1),
	    fdGetPixelCharge(fd, x+1, y) + fdGetPixelCharge(fd, x, y-1)	    
	    };*/
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


/* Deprecated
void fdMarkEvents(eROSITAEvent* list, int* nlist, 
		  int* maxidx, int* minidx,
		  FramestoreDetector* fd, 
		  int x, int y)
{
  // TODO Determine split threshold from charge of MAXIMUM split partner.

  // Split threshold, depending on the charge in the main pixel.
  const double split_threshold = 
    fd->pixels.array[x][y].charge * fd->split_threshold;

  printf("split threshold: %lf\n", split_threshold);

  // Possible neighbors (no diagonal neighbors).
  const int neighbors_x[4] = {0, 0, 1, -1};
  const int neighbors_y[4] = {1, -1, 0, 0};

  int neighbor, xi, yi;

  // Create a new event in the list.
  assert(*nlist+1 < MAX_N_SPLIT_LIST);  
  list[*nlist].pha = getChannel(fd->pixels.array[x][y].charge, 
				fd->generic.rmf);
  list[*nlist].energy = fd->pixels.array[x][y].charge * 1.e3; // [eV]
  list[*nlist].xi  = x;
  list[*nlist].yi  = y;
  (*nlist)++;
  // Delete the charge of this event from the pixel array.
  fd->pixels.array[x][y].charge = 0.;

  // Check if the new event has the maximum or minium energy in the 
  // event list.
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
*/


void fdPatternIdentification(eROSITAEvent* list, const int nlist, 
			     const int maxidx, const int minidx)
{
 
  // Single events.
  if (1==nlist) {
    list[0].pat_inf = 0;
  }

  // Double events.
  else if (2==nlist) {
    // 0i0
    // 0x0
    // 000
    if ((list[maxidx].xi==list[minidx].xi) && 
	(list[maxidx].yi==list[minidx].yi-1)) {
      list[maxidx].pat_inf = 15;
      list[minidx].pat_inf = 12;
    }
    // 000
    // 0xi
    // 000
    else if ((list[maxidx].xi==list[minidx].xi-1) && 
	     (list[maxidx].yi==list[minidx].yi)) {
      list[maxidx].pat_inf = 25;
      list[minidx].pat_inf = 26;
    }
    // 000
    // 0x0
    // 0i0
    else if ((list[maxidx].xi==list[minidx].xi) && 
	     (list[maxidx].yi==list[minidx].yi+1)) {
      list[maxidx].pat_inf = 35;
      list[minidx].pat_inf = 38;
    }
    // 000
    // ix0
    // 000
    else if ((list[maxidx].xi==list[minidx].xi+1) && 
	     (list[maxidx].yi==list[minidx].yi)) {
      list[maxidx].pat_inf = 45;
      list[minidx].pat_inf = 44;
    }
    // Invalid pattern. (Actually this code should never be executed.)
    else {
      headas_chat(5, "Invalid double event: (%d,%d) and (%d,%d)!\n",
		  list[maxidx].xi, list[maxidx].yi,
		  list[minidx].xi, list[minidx].yi);
      list[maxidx].pat_inf = FD_INVALID_PATTERN;
      list[minidx].pat_inf = FD_INVALID_PATTERN;
    }
  } // END of double events.
  
  // Triple events.
  else if (3==nlist) {

    // Find the index of the non-maximum and non-minimum split partner.
    int mididx;
    for (mididx=0; mididx<nlist; mididx++) {
      if ((mididx!=maxidx) && (mididx!=minidx)) break;
    }

    // Check if it is a valid triple split pattern, i.e., around the corner
    // with the maximum event in the the corner.
    if ((abs(list[maxidx].xi-list[minidx].xi) + abs(list[maxidx].xi-list[mididx].xi) == 1) &&
	(abs(list[maxidx].yi-list[minidx].yi) + abs(list[maxidx].yi-list[mididx].yi) == 1)) {

      // Check the orientation of the triple split pattern.
      // 0i0
      // 0xm
      // 000
      if ((list[maxidx].xi==list[mididx].xi-1) && 
	  (list[maxidx].yi==list[minidx].yi-1)) {
	list[maxidx].pat_inf = 55;
	list[mididx].pat_inf = 56;
	list[minidx].pat_inf = 52;      
      }

      // 0m0
      // 0xi
      // 000
      else if ((list[maxidx].xi==list[minidx].xi-1) && 
	       (list[maxidx].yi==list[mididx].yi-1)) {
	list[maxidx].pat_inf = 55;
	list[mididx].pat_inf = 52;
	list[minidx].pat_inf = 56;      
      }

      // 000
      // 0xi
      // 0m0
      else if ((list[maxidx].xi==list[minidx].xi-1) && 
	       (list[maxidx].yi==list[mididx].yi+1)) {
	list[maxidx].pat_inf = 65;
	list[mididx].pat_inf = 68;
	list[minidx].pat_inf = 66;      
      }
     
      // 000
      // 0xm
      // 0i0
      else if ((list[maxidx].xi==list[mididx].xi-1) && 
	       (list[maxidx].yi==list[minidx].yi+1)) {
	list[maxidx].pat_inf = 65;
	list[mididx].pat_inf = 66;
	list[minidx].pat_inf = 68;      
      }

      // 000
      // mx0
      // 0i0
      else if ((list[maxidx].xi==list[mididx].xi+1) && 
	       (list[maxidx].yi==list[minidx].yi+1)) {
	list[maxidx].pat_inf = 75;
	list[mididx].pat_inf = 74;
	list[minidx].pat_inf = 78;      
      }

      // 000
      // ix0
      // 0m0
      else if ((list[maxidx].xi==list[minidx].xi+1) && 
	       (list[maxidx].yi==list[mididx].yi+1)) {
	list[maxidx].pat_inf = 75;
	list[mididx].pat_inf = 78;
	list[minidx].pat_inf = 74;      
      }

      // 0m0
      // ix0
      // 000
      else if ((list[maxidx].xi==list[minidx].xi+1) && 
	       (list[maxidx].yi==list[mididx].yi-1)) {
	list[maxidx].pat_inf = 85;
	list[mididx].pat_inf = 82;
	list[minidx].pat_inf = 84;      
      }

      // 0i0
      // mx0
      // 000
      else if ((list[maxidx].xi==list[mididx].xi+1) && 
	       (list[maxidx].yi==list[minidx].yi-1)) {
	list[maxidx].pat_inf = 85;
	list[mididx].pat_inf = 84;
	list[minidx].pat_inf = 82;      
      }

      // Invalid pattern. (Actually this code should never be executed.)
      else {
	headas_printf("Invalid triple event!\n");
	list[maxidx].pat_inf = FD_INVALID_PATTERN;
	list[mididx].pat_inf = FD_INVALID_PATTERN;
	list[minidx].pat_inf = FD_INVALID_PATTERN;
      }

    } else { // It is not valid triple pattern with the maximum in the corner.
      list[maxidx].pat_inf = FD_INVALID_PATTERN;
      list[mididx].pat_inf = FD_INVALID_PATTERN;
      list[minidx].pat_inf = FD_INVALID_PATTERN;
      headas_chat(5, "Invalid triple event: max (%d,%d)  mid (%d,%d) min (%d,%d)\n",
		  list[maxidx].xi, list[maxidx].yi, 
		  list[mididx].xi, list[mididx].yi,
		  list[minidx].xi, list[minidx].yi);
    }
  } // END of triple events.

  // Quadruple events.
  else if (4==nlist) {
    int count;

    // Find the 2 events, which are neither the maximum nor the minimum.
    int mididx[2] = { maxidx, maxidx };
    for (count=0; count<nlist; count++) {
      if ((count!=maxidx) && (count!=minidx)) {
	if (mididx[0] == maxidx) {
	  mididx[0] = count;
	} else {
	  mididx[1] = count;
	}
      }
    }

    // Check if the 4 events form a valid 2x2 matrix with the maximum
    // and minimum events on opposite corners.
    int check=0;
    for (count=0; count<2; count++) {
      if (((list[mididx[count]].xi == list[maxidx].xi) &&
	   (list[mididx[count]].yi == list[minidx].yi)) ||
	  ((list[mididx[count]].yi == list[maxidx].yi) &&
	   (list[mididx[count]].xi == list[minidx].xi))) {
	check++;
      }
    }

    if (2!=check) {
      // Invalid pattern.
      for (count=0; count<nlist; count++) {
	list[count].pat_inf = FD_INVALID_PATTERN;
      }
      return;
    } 

    // The pattern is a valid 2x2 matrix.
    // Check the orientation of the pattern.
    // 0mi
    // 0xm
    // 000
    if ((list[maxidx].xi == list[minidx].xi-1) &&
	(list[maxidx].yi == list[minidx].yi-1)) {
      list[maxidx].pat_inf = 95;
      list[minidx].pat_inf = 93;
      
      if (list[maxidx].xi == list[mididx[0]].xi) {
	list[mididx[0]].pat_inf = 92;
	list[mididx[1]].pat_inf = 96;
      } else {
	list[mididx[0]].pat_inf = 96;
	list[mididx[1]].pat_inf = 92;
      }
    }

    // 000
    // 0xm
    // 0mi
    else if ((list[maxidx].xi == list[minidx].xi-1) &&
	     (list[maxidx].yi == list[minidx].yi+1)) {
      list[maxidx].pat_inf = 105;
      list[minidx].pat_inf = 109;
      
      if (list[maxidx].xi == list[mididx[0]].xi) {
	list[mididx[0]].pat_inf = 108;
	list[mididx[1]].pat_inf = 106;
      } else {
	list[mididx[0]].pat_inf = 106;
	list[mididx[1]].pat_inf = 108;
      }
    }

    // 000
    // mx0
    // im0
    else if ((list[maxidx].xi == list[minidx].xi+1) &&
	     (list[maxidx].yi == list[minidx].yi+1)) {
      list[maxidx].pat_inf = 115;
      list[minidx].pat_inf = 117;
      
      if (list[maxidx].xi == list[mididx[0]].xi) {
	list[mididx[0]].pat_inf = 118;
	list[mididx[1]].pat_inf = 114;
      } else {
	list[mididx[0]].pat_inf = 114;
	list[mididx[1]].pat_inf = 118;
      }
    }

    // im0
    // mx0
    // 000
    else if ((list[maxidx].xi == list[minidx].xi+1) &&
	     (list[maxidx].yi == list[minidx].yi-1)) {
      list[maxidx].pat_inf = 125;
      list[minidx].pat_inf = 121;
      
      if (list[maxidx].xi == list[mididx[0]].xi) {
	list[mididx[0]].pat_inf = 122;
	list[mididx[1]].pat_inf = 124;
      } else {
	list[mididx[0]].pat_inf = 124;
	list[mididx[1]].pat_inf = 122;
      }
    }

    // Invalid pattern. (Actually this code should never be executed.)
    else {
      headas_printf("Invalid quadruple event!\n");
      for (count=0; count<nlist; count++) {
	list[count].pat_inf = FD_INVALID_PATTERN;
      }
    }

  } // END of quadruple events.

  // Invalid patterns with more than 4 events.
  else {
    int count;
    for(count=0; count<nlist; count++) {
      list[count].pat_inf = FD_INVALID_PATTERN;
    }
  } // END of invalid patterns with more than 4 events.

}

