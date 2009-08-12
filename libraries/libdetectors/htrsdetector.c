#include "htrsdetector.h"


int initHTRSDetector(HTRSDetector* hd, struct HTRSDetectorParameters* parameters)
{
  int status = EXIT_SUCCESS;

  // Call the initialization routines of the underlying data structures.
  status = initGenericDetector(&hd->generic, &parameters->generic);
  if (EXIT_SUCCESS!=status) return(status);
  status = initHexagonalPixels(&hd->pixels, &parameters->pixels);
  if (EXIT_SUCCESS!=status) return(status);

  // Set up the HTRS configuration.
  hd->dead_time = parameters->dead_time;

  // Create a new FITS event file and open it.
  status = openNewHTRSEventFile(&hd->eventlist, parameters->eventlist_filename,
				parameters->eventlist_template);
  if (EXIT_SUCCESS!=status) return(status);

  return(status);
}



int cleanupHTRSDetector(HTRSDetector* hd)
{
  int status=EXIT_SUCCESS;

  // Call the cleanup routines of the underlying data structures.
  cleanupHexagonalPixels(&hd->pixels);
  status = closeHTRSEventFile(&hd->eventlist);

  return(status);
}



int addImpact2HTRSDetector(HTRSDetector* hd, Impact* impact)
{
  int status=EXIT_SUCCESS;

  // Determine a detector channel (PHA channel) according to the RMF.
  // The channel is obtained from the RMF using the corresponding
  // HEAdas routine which is based on drawing a random number.
  long channel;
  ReturnChannel(hd->generic.rmf, impact->energy, 1, &channel);

  // Check if the photon is really measured. If the
  // PHA channel returned by the HEAdas RMF function is equal to '-1', 
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
  float charge = getEnergy(channel, hd->generic.rmf);
  
  if (charge > 0.) {

    // Determine the affected pixel(s) and their corresponding charge fractions.
    int pixel[2];
    double fraction[2];
    int nsplits = getHexagonalPixelSplits(&hd->pixels, &hd->generic, impact->position, 
					  pixel, fraction);

    // Loop over all split partners.
    int pixel_counter;
    for (pixel_counter=0; pixel_counter<nsplits; pixel_counter++) {
      if (INVALID_PIXEL != pixel[pixel_counter]) {
	
	// Check if the pixel is currently active.
	if (impact->time - hd->pixels.array[pixel[pixel_counter]].last_impact < hd->dead_time) {
	  // The photon cannot be detected in this pixel as it is still within the
	  // dead time after the previous event.
	  continue;
	}

	// Add the charge created by the photon to the affected detector pixel(s).
	HTRSEvent event;

	// Determine the detector channel that corresponds to the charge fraction
	// created by the incident photon in the regarded pixel.
	event.pha = getChannel(charge * fraction[pixel_counter], hd->generic.rmf);
	//                     |        |-> charge fraction due to split events
	//                     |-> charge created by incident photon
      
	// Check lower threshold (PHA and energy):
	if ((event.pha>=hd->generic.pha_threshold) && 
	    (charge*fraction[pixel_counter]>=hd->generic.energy_threshold)) { 
	
	  // TODO REMOVE
	  assert(event.pha >= 0);
	  // Maybe: if (event.pha < 0) continue;
	  
	  // The impact has generated an event in this pixel,
	  // so add it to the event file.
	  event.pixel = pixel[pixel_counter];
	  event.time = impact->time;
	  
	  // Store the time of this impact in the pixel in order to consider the
	  // dead time.
	  hd->pixels.array[pixel[pixel_counter]].last_impact = impact->time;

	  // Add event to event file.
	  status = addHTRSEvent2File(&hd->eventlist, &event);
	  if (EXIT_SUCCESS!=status) return(status);

	} // END Check for thresholds.
      } // END if(charge>0.)
    } // END if valid pixel
  } // END of loop over split partners

  return(status);
}


