#include "wfidetector.h"


////////////////////////////////////////////////////////////////////
// Static function declarations
////////////////////////////////////////////////////////////////////


/** Check the detector sensitivity in the given detector line at the
    specified time. This routine checks, whether the WFIDetector is
    sensitive in line 'y' at the given time. According to the
    currently implemented readout strategy the WFIDetector is
    insensitive to new photons in an interval immediately after the
    line readout. This interval is called 'line_clear_time'. This
    function checks whether the designated line was just read out and
    the time lies within the clear time.  It returns 1 if the detector
    is sensitive in the line 'y' at the specified time. Otherwise the
    return value is '0'. */
static inline int WFIDetectorIsSensitive(int y, WFIDetector*, double time);

/** Read out the charge from the currently active readout lines (one
    or two) of the WFI detector. This function implements the
    different readout modes of the WFIDetector model. According to
    the chosen number of readout lines it performs the readout of the
    currently active detector lines. The events read from the
    detector pixels are stored in the output even file. */
static inline int readoutLinesWFIDetector(WFIDetector*);


////////////////////////////////////////////////////////////////////
// Program Code
////////////////////////////////////////////////////////////////////


int initWFIDetector(WFIDetector* wd, struct WFIDetectorParameters* parameters)
{
  int status = EXIT_SUCCESS;

  // Call the initialization routines of the underlying data structures.
  status = initGenericDetector(&wd->generic, &parameters->generic);
  if (EXIT_SUCCESS!=status) return(status);
  wd->pixels = newSquarePixels(&parameters->pixels, &status);
  if (EXIT_SUCCESS!=status) return(status);

  // Set up the WFI configuration.
  wd->line_readout_time = parameters->line_readout_time;
  wd->line_clear_time   = parameters->line_clear_time;
  wd->readout_directions = parameters->readout_directions;

  // Set the first readout time such that the first readout is performed 
  // immediately at the beginning of the simulation.
  wd->readout_time = parameters->t0 - wd->line_readout_time;
  wd->readout_lines[0] = 0;
  wd->readout_lines[1] = wd->pixels->ywidth-1;
  wd->frame = -1;

  // Create a new event list FITS file and open it.
  status = openNewWFIEventFile(&wd->eventlist, parameters->eventlist_filename,
				   parameters->eventlist_template);
  if (EXIT_SUCCESS!=status) return(status);

  return(status);
}



int cleanupWFIDetector(WFIDetector* wd)
{
  int status=EXIT_SUCCESS;

  // Call the cleanup routines of the underlying data structures.
  destroySquarePixels(&wd->pixels);
  cleanupGenericDetector(&wd->generic);
  status = closeWFIEventFile(&wd->eventlist);

  return(status);
}



int addImpact2WFIDetector(WFIDetector* wd, Impact* impact)
{
  int status=EXIT_SUCCESS;

  // Before adding the new impact to the detector check whether
  // a readout has to be performed in advance.
  status=checkReadoutWFIDetector(wd, impact->time);

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
    return(status); // Break the function (and continue with the next photon).
  }
  assert(channel>=0);

  // Get the corresponding created charge.
  // NOTE: In this simulation the charge is represented by the nominal
  // photon energy which corresponds to the PHA channel according to the
  // EBOUNDS table.
  float charge = getEBOUNDSEnergy(channel, wd->generic.rmf, 0);
  
  if (charge > 0.) {
    int x[4], y[4];
    double fraction[4];
    
    // Determine the affected detector pixels.
    int npixels = 
      getSquarePixelsGaussianSplits(wd->pixels, &(wd->generic.gcc), 
				    impact->position, x, y, fraction);
    
    // Set the valid flag for the affected pixel in order to 
    // keep a record, whether the pattern has been generated
    // by a single photon.
    SPupdateValidFlag(wd->pixels, x, y, npixels);
  
    // Add the charge created by the photon to the affected detector pixels.
    int count;
    for (count=0; count<npixels; count++) {
      SPaddCharge(wd->pixels, x[count], y[count], 
		  charge * fraction[count] * 
		  // |      |-> charge fraction due to split events
		  // |-> charge created by incident photon
		  WFIDetectorIsSensitive(y[count], wd, impact->time));
		  // |-> "1" if pixel can measure charge, "0" else
    }
  } // END if(charge>0.)

  return(status);
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
	(1==wd->readout_directions)?(wd->pixels->ywidth-1):(wd->pixels->yoffset-1);
      wd->readout_lines[1] = 
	(1==wd->readout_directions)?(-1):(wd->pixels->yoffset);

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
      if (wd->pixels->line2readout[wd->readout_lines[lineindex]]>0) {    
	clearLineSquarePixels(wd->pixels, wd->readout_lines[lineindex]);
      }
    }
  }
  
  return(status);
}



static inline int readoutLinesWFIDetector(WFIDetector* wd)
{
  int x, lineindex;
  int status=EXIT_SUCCESS;

  // Read out only the 1 or 2 detector read-out lines from the whole detector array.
  for (lineindex=0; lineindex<wd->readout_directions; lineindex++) {
    // Check whether any pixel in this line has been affected by an event since
    // the last read-out cycle. If not, the line can be neglected.
    if (1==wd->pixels->line2readout[wd->readout_lines[lineindex]]) {
      for (x=0; x<wd->pixels->xwidth; x++) {
	if (wd->pixels->array[x][wd->readout_lines[lineindex]].charge > 1.e-6) {
	  // Determine the detector channel that corresponds to the charge stored
	  // in the detector pixel.
	  WFIEvent event = {
	    .pha = getEBOUNDSChannel(wd->pixels->array[x][wd->readout_lines[lineindex]].charge, 
				     wd->generic.rmf)
	  };

	  // Check lower threshold (PHA and energy):
	  if ((event.pha>=wd->generic.pha_threshold) && 
	      (wd->pixels->array[x][wd->readout_lines[lineindex]].charge>=
	       wd->generic.energy_threshold)) { 
	    
	    // TODO REMOVE
	    assert(event.pha >= 0);
	    // Maybe: if (event.pha < 0) continue;
	    
	    // There is an event in this pixel, so insert it into the eventlist:
	    event.time = wd->readout_time;
	    event.xi = x;
	    event.yi = wd->readout_lines[lineindex];
	    event.frame = wd->frame;
	    event.patnum = 0;
	    event.patid = -1;
	    event.pileup = 0;
	    event.f_valid = 
	      wd->pixels->array[x][wd->readout_lines[lineindex]].valid_flag;
	    
	    status = addWFIEvent2File(&wd->eventlist, &event);
	    if (EXIT_SUCCESS!=status) return(status);
	  } // END of check for threshold
	} // END of check whether  charge > 1.e-6
      } // END of loop over x
    } // END of check whether any pixel in the line has been affected by an impact
  } // END of loop over readout lines

  return(status);
}



static inline int WFIDetectorIsSensitive(int y, WFIDetector* wd, double time)
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



