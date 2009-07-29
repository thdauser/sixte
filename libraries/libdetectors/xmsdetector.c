#include "xmsdetector.h"


int initXMSDetector(XMSDetector* xd, struct XMSDetectorParameters* parameters)
{
  int status = EXIT_SUCCESS;

  // Call the initialization routines of the underlying data structures.
  status = initGenericDetector(&xd->generic, &parameters->generic);
  if (EXIT_SUCCESS!=status) return(status);
  //  status = initSquarePixels(&xd->pixels, &parameters->pixels);
  if (EXIT_SUCCESS!=status) return(status);

  // Set up the XMS configuration.
  // --- Currently nothing to do. ---

  // Create a new FITS event file and open it.
  status = openNewXMSEventFile(&xd->eventlist, parameters->eventlist_filename,
			       parameters->eventlist_template);
  if (EXIT_SUCCESS!=status) return(status);

  return(status);
}



int cleanupXMSDetector(XMSDetector* xd)
{
  int status=EXIT_SUCCESS;

  // Call the cleanup routines of the underlying data structures.
  cleanupSquarePixels(&xd->pixels);
  status = closeXMSEventFile(&xd->eventlist);

  return(status);
}




int addImpact2XMSDetector(XMSDetector* xd, Impact* impact)
{
  int status=EXIT_SUCCESS;

  // Determine a detector channel (PHA channel) according to the RMF.
  // The channel is obtained from the RMF using the corresponding
  // HEAdas routine which is based on drawing a random number.
  long channel;
  ReturnChannel(xd->generic.rmf, impact->energy, 1, &channel);

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
  float charge = getEnergy(channel, xd->generic.rmf);
  
  if (charge > 0.) {
    int x[4], y[4];
    double fraction[4];
    
    // Determine the affected detector pixels.
    int npixels = getSquarePixelsSplits(&xd->pixels, &xd->generic, impact->position, 
					x, y, fraction);
    
    // Add the charge created by the photon to the affected detector pixels.
    int count;
    for (count=0; count<npixels; count++) {
      if (x[count] != INVALID_PIXEL) {
	XMSEvent event;

	// Determine the detector channel that corresponds to the charge fraction
	// created by the incident photon in the regarded pixel.
	event.pha = getChannel(charge * fraction[count], xd->generic.rmf);
	//                     |        |-> charge fraction due to split events
	//                     |-> charge created by incident photon

	// Check lower threshold (PHA and energy):
	if ((event.pha>=xd->generic.pha_threshold) && 
	    (charge*fraction[count]>=xd->generic.energy_threshold)) { 

	  // TODO REMOVE
	  assert(event.pha >= 0);
	  // Maybe: if (event.pha < 0) continue;

	  // The impact has generated an event in this pixel,
	  // so add it to the event file.
	  event.xi = x[count];
	  event.yi = y[count];
	  event.time = impact->time;

	  // Add event to event file.
	  status = addXMSEvent2File(&xd->eventlist, &event);
	  if (EXIT_SUCCESS!=status) return(status);

	} // END Check for thresholds.
      } // END if valid pixel
    } // END of loop over all affected pixels.
  } // END if(charge>0.)

  return(status);
}


