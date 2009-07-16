#include "wfidetector.h"


int initWFIDetector(WFIDetector* wd, struct WFIDetectorParameters* parameters)
{
  int status = EXIT_SUCCESS;

  // Call the initialization routines of the underlying data structures.
  status = initGenericDetector(&wd->generic, &parameters->generic);
  if (EXIT_SUCCESS!=status) return(status);
  status = initSquarePixels(&wd->pixels, &parameters->pixels);
  if (EXIT_SUCCESS!=status) return(status);

  // Set up the WFI configuration.
  wd->line_readout_time = parameters->line_readout_time;
  wd->line_clear_time   = parameters->line_clear_time;
  wd->readout_directions = parameters->readout_directions;

  // Set the first readout time such that the first readout is performed 
  // immediately at the beginning of the simulation.
  wd->readout_time = parameters->t0 - wd->line_readout_time;
  wd->readout_lines[0] = 0;
  wd->readout_lines[1] = wd->pixels.ywidth-1;
  wd->frame = -1;

  // Create a new event list FITS file and open it.
  status = openNewWFIEventFile(&wd->eventlist, parameters->eventlist_filename,
				   parameters->eventlist_template);
  if (EXIT_SUCCESS!=status) return(status);

  return(status);
}




void cleanupWFIDetector(WFIDetector* wd)
{
  // TODO: Perform one last readout !!

  // Call the cleanup routines of the underlying data structures.
  cleanupSquarePixels(&wd->pixels);
  closeWFIEventFile(&wd->eventlist);
}



void addImpact2WFIDetector(WFIDetector* wd, Impact* impact)
{
  // Determine a detector channel (PHA channel) according to the RMF.
  // The channel is obtained from the RMF using the corresponding
  // HEAdas routine which is based on drawing a random number.
  long channel;
  ReturnChannel(wd->generic.rmf, impact->energy, 1, &channel);

  // Check if the photon is really measured. If the
  // PHA channel returned by the HEAdas RMF function is '-1', 
  // the photon is not detected.
  // This can happen as the RMF is actually an RSP including, e.g., 
  // the detector quantum efficiency and filter transmission.
  if (-1==channel) {
    return; // Break the function (and continue with the next photon).
  }
  assert(channel>=0);

  // Get the corresponding created charge.
  // NOTE: In this simulation the charge is represented by the nominal
  // photon energy which corresponds to the PHA channel according to the
  // EBOUNDS table.
  float charge = getEnergy(channel, wd->generic.rmf);
  
  if (charge > 0.) {
    int x[4], y[4];
    double fraction[4];
    
    // Determine the affected detector pixels.
    int npixels = getSquarePixelsSplits(&wd->pixels, &wd->generic, impact->position, 
					x, y, fraction);
    
    // Add the charge created by the photon to the affected detector pixels.
    int count;
    for (count=0; count<npixels; count++) {
      if (x[count] != INVALID_PIXEL) {
	wd->pixels.array[x[count]][y[count]].charge += 
	  charge * fraction[count] * 
	  // |      |-> charge fraction due to split events
	  // |-> charge created by incident photon
	  WFIDetectorIsSensitive(x[count], y[count], wd, impact->time);
	  // |-> "1" if pixel can measure charge, "0" else
      }
    }
  } // END if(charge>0.)
}




int checkReadoutWFIDetector(WFIDetector* wd, double time) 
{
  int status=EXIT_SUCCESS;

  // The WFI detector is read out line by line, with two readout lines 
  // starting in the middle of the detector array. The readout of one 
  // individual line requires the 'line_readout_time'. The readout is performed at 
  // the beginning of this interval. Then the charges in the line are cleared. 
  // During the 'line_clear_time' the pixels in the detector line are inactive, i.e., 
  // they cannot receive new charges from incident photons during that time.
  // Photons that hit the pixel during the readout time but after the clear time are 
  // accepted and the created charge is stored in the pixel until the next readout
  // process.

  // Determine, which line / which 2 lines currently have to be read out.
  while (time > wd->readout_time + wd->line_readout_time) {
    // The current line number has to be decreased (readout process starts in the
    // middle of the detector).
    wd->readout_lines[0]--;
    wd->readout_lines[1]++;
    if (wd->readout_lines[0] < 0) {       
      // Start new detector frame:
      wd->frame++; 
      wd->readout_lines[0] = 
	(1==wd->readout_directions)?(wd->pixels.ywidth-1):(wd->pixels.yoffset-1);
      wd->readout_lines[1] = 
	(1==wd->readout_directions)?(-1):(wd->pixels.yoffset);

      // Print the time of the current events in order (program status
      // information for the user).
      if (0 == wd->frame % 100000) {
	// Display this status information only each 100 000 frames in order to 
	// avoid numerical load due to the displaying on STDOUT.
	headas_printf("\rtime: %.3lf s ", wd->readout_time);
	fflush(NULL);
      }
    }

    // Update the current detector readout time, which is used in the 
    // event list output.
    wd->readout_time += wd->line_readout_time;

    // Perform the readout on the 1 or 2 (!) current lines 
    // (i.e., the one or two new lines, chosen in the
    // step before) and write the data to the FITS file.
    // After the readout clear the  lines.

    status = readoutLinesWFIDetector(wd);
    if(EXIT_SUCCESS!=status) return(status);

    // Clear the previously read out lines.
    int lineindex;
    for(lineindex=0; lineindex<wd->readout_directions; lineindex++) {
      clearLineSquarePixels(&wd->pixels, wd->readout_lines[lineindex]);
    }
  }
  
  return(status);
}





inline int readoutLinesWFIDetector(WFIDetector* wd)
{
  int x, lineindex;
  int status = EXIT_SUCCESS;

  // Read out the entire detector array.
  for (lineindex=0; lineindex<wd->readout_directions; lineindex++) {
    for (x=0; x<wd->pixels.xwidth; x++) {
      if (wd->pixels.array[x][wd->readout_lines[lineindex]].charge > 1.e-6) {
	WFIEvent event;
	// Determine the detector channel that corresponds to the charge stored
	// in the detector pixel.
	event.pha = getChannel(wd->pixels.array[x][wd->readout_lines[lineindex]].charge, 
			       wd->generic.rmf);

	// Check lower threshold (PHA and energy):
	if ((event.pha>=wd->generic.pha_threshold) && 
	    (wd->pixels.array[x][wd->readout_lines[lineindex]].charge>=
	     wd->generic.energy_threshold)) { 

	  // REMOVE TODO
	  assert(event.pha >= 0);
	  // Maybe: if (event.pha < 0) continue;
	
	  // There is an event in this pixel, so insert it into the eventlist:
	  event.time = wd->readout_time;
	  event.xi = x;
	  event.yi = wd->readout_lines[lineindex];
	  event.frame = wd->frame;
	  status = addWFIEvent2File(&wd->eventlist, &event);
	  if (EXIT_SUCCESS!=status) return(status);
	} // END of check for threshold
      } // END of check whether  charge > 1.e-6
    } // END of loop over x
  } // END of loop over readout lines

  return(status);
}




inline int WFIDetectorIsSensitive(int x, int y, WFIDetector* wd, double time)
{
  if ((y==wd->readout_lines[0])||
      ((2==wd->readout_directions)&&(y==wd->readout_lines[1]))) {
    // If we are at the beginning of the readout interval of the regarded line
    // within the clear time, the photon is not measured.
    if (time - wd->readout_time < wd->line_clear_time) {
      // The specified line is cleared at the moment, so no photon can 
      // be detected.
      return(0);
    }
  }
  
  // The specified detector pixel is active at the moment.
  return(1);
}



