#include "xmsdetector.h"


int initXMSDetector(XMSDetector* xd, struct XMSDetectorParameters* parameters)
{
  int status = EXIT_SUCCESS;

  // Call the initialization routines of the underlying data structures.
  // Inner TES pixel array:
  status = initGenericDetector(&xd->generic_inner, &parameters->generic_inner);
  if (EXIT_SUCCESS!=status) return(status);
  status = initSquarePixels(&xd->pixels_inner, &parameters->pixels_inner);
  if (EXIT_SUCCESS!=status) return(status);

  // Outer TES pixel array:
  status = initGenericDetector(&xd->generic_outer, &parameters->generic_outer);
  if (EXIT_SUCCESS!=status) return(status);
  status = initSquarePixels(&xd->pixels_outer, &parameters->pixels_outer);
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
  cleanupSquarePixels(&xd->pixels_inner);
  cleanupSquarePixels(&xd->pixels_outer);
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
  ReturnChannel(xd->generic_inner.rmf, impact->energy, 1, &channel);

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
  float charge = getEnergy(channel, xd->generic_inner.rmf);
  
  if (charge > 0.) {
    int x[4], y[4];
    double fraction[4];

    // Check whether the impact position lies within the inner or the outer TES
    // microcalorimeter array.
    XMSEvent event={ .grade=0 };
    int npixels;
    if ((impact->position.x > -xd->pixels_inner.xoffset*xd->pixels_inner.xpixelwidth) &&
	(impact->position.x <  xd->pixels_inner.xoffset*xd->pixels_inner.xpixelwidth) &&
	(impact->position.y > -xd->pixels_inner.yoffset*xd->pixels_inner.ypixelwidth) &&
	(impact->position.y <  xd->pixels_inner.yoffset*xd->pixels_inner.ypixelwidth)) {
      // Impact position is within the inner pixel array.
      event.array = 1;
      // Determine the affected detector pixels.
      npixels = getSquarePixelsGaussianSplits(&xd->pixels_inner, &(xd->generic_inner.gcc), 
					      impact->position, x, y, fraction);
    } else { 
      // Impact position is either within the outer pixel array or even completely
      // outside the detector.
      event.array = 2;
      // Determine the affected detector pixels.
      npixels = getSquarePixelsGaussianSplits(&xd->pixels_outer, &(xd->generic_outer.gcc), 
					      impact->position, x, y, fraction);
    }

    
    // Add the charge created by the photon to the affected detector pixels.
    int count;
    for (count=0; count<npixels; count++) {
      if (x[count] != INVALID_PIXEL) {

	long pha_threshold;
	float energy_threshold;
	// Determine the detector channel that corresponds to the charge fraction
	// created by the incident photon in the regarded pixel.
	if (1 == event.array) {
	  event.pha = getChannel(charge * fraction[count], xd->generic_inner.rmf);
	  //                     |        |-> charge fraction due to split events
	  //                     |-> charge created by incident photon
	  pha_threshold = xd->generic_inner.pha_threshold;
	  energy_threshold = xd->generic_inner.energy_threshold;
	} else {
	  event.pha = getChannel(charge * fraction[count], xd->generic_outer.rmf);
	  //                     |        |-> charge fraction due to split events
	  //                     |-> charge created by incident photon
	  pha_threshold = xd->generic_outer.pha_threshold;
	  energy_threshold = xd->generic_outer.energy_threshold;
	}

	// Check lower threshold (PHA and energy):
	if ((event.pha>=pha_threshold) && (charge*fraction[count]>=energy_threshold)) { 

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
      }// END if valid pixel
    } // END of loop over all affected pixels.
  } // END if(charge>0.)

  return(status);
}


