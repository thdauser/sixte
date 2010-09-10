#include "wfi_pattern_recombination.h"


int pattern_recombination_getpar(struct Parameters* parameters)
{
  int status = EXIT_SUCCESS;

  if ((status = PILGetFname("eventlist_filename", parameters->eventlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the input file!\n", status);
  }

  else if ((status = PILGetFname("pattern_filename", parameters->pattern_filename))) {
    HD_ERROR_THROW("Error reading the name of the output file!\n", status);
  }

  else if ((status = PILGetFname("response_filename", parameters->response_filename))) {
    HD_ERROR_THROW("Error reading the name of the detector response file!\n", status);
  }

  // Get the name of the FITS template directory.
  // First try to read it from the environment variable.
  // If the variable does not exist, read it from the PIL.
  else { 
    char* buffer;
    if (NULL!=(buffer=getenv("SIXT_FITS_TEMPLATES"))) {
      strcpy(parameters->templatedir, buffer);
    } else {
      if ((status = PILGetFname("templatedir", parameters->templatedir))) {
	HD_ERROR_THROW("Error reading the path of the FITS templates!\n", status);
      }
    }
  }

  return(status);
}



WFIEvent mark(int row, int column, int rows, int columns, WFIEvent** ccdarr, 
	      WFIEvent* maxevent, WFIEvent* evtlist, long* nevtlist, struct RMF* rmf) 
{
  // Return immediately if the pixel is empty.
  if (-1 == ccdarr[column][row].patnum) {
    return(ccdarr[column][row]);
  }

  // Remember this event ...
  WFIEvent myevent = ccdarr[column][row];
  myevent.patnum = 1;
  
  // ... and mark as deleted.
  ccdarr[column][row].patnum = -1;

  evtlist[(*nevtlist)++] = myevent;
  if (myevent.pha > maxevent->pha) {
    *maxevent = myevent;
  }

  // Visit neighbours:
  int visitrow[8]={row-1 , row+1 , row     , row     , row-1   , row-1   , row+1   , row+1   };
  int visitcol[8]={column, column, column-1, column+1, column-1, column+1, column-1, column+1};

  int i;
  for (i=0; i<8; i++) {
    // Border event?
    if ((visitrow[i] < 0) || (visitrow[i] >= rows) ||
	(visitcol[i] < 0) || (visitcol[i] >= columns)) {
      myevent.patnum = 1000;
    } else {
      // Normal events:
      WFIEvent event = mark(visitrow[i], visitcol[i], rows, columns, ccdarr,
			    maxevent, evtlist, nevtlist, rmf);
      if (-1 != event.patnum ) {
	myevent.pileup += event.pileup; // Remember if pileup occurred in subevent
	myevent.patnum += event.patnum;
	myevent.pha = getChannel(getEnergy(myevent.pha, rmf, 0)+
				 getEnergy(event.pha, rmf, 0), rmf);
      }
    }
  } // End of loop over all neighbours

  // Return combined event:
  myevent.xi = maxevent->xi;
  myevent.yi = maxevent->yi;

  return(myevent);

} // End of mark()



long min(long* array, long nelements)
{
  long minidx=0;
  long i;
  for (i=0; i<nelements; i++) {
    if (array[minidx] > array[i]) {
      minidx = i;
    }
  }
  return(array[minidx]);
}



long max(long* array, long nelements)
{
  long maxidx=0;
  long i;
  for (i=0; i<nelements; i++) {
    if (array[maxidx] < array[i]) {
      maxidx = i;
    }
  }
  return(array[maxidx]);
}



WFIEvent pattern_id(WFIEvent* components, long ncomponents, WFIEvent event)
{
  long i;

  // Get events with maximum and minimum energy:
  long maximum=0, minimum=1000000;
  long nmaxidx=0, nminidx=0;
  long amaxidx[ARRAY_LENGTH], aminidx[ARRAY_LENGTH];
  // Initialize with 0 values.
  for(nmaxidx=0; nmaxidx<ARRAY_LENGTH; nmaxidx++) {
    amaxidx[nmaxidx]=0;
    aminidx[nmaxidx]=0;
  }
  nmaxidx=0;
  for (i=0; i<ncomponents; i++) {
    if (components[i].pha > maximum) maximum = components[i].pha;
    if (components[i].pha < minimum) minimum = components[i].pha;
  }
  for (i=0; i<ncomponents; i++) {
    if (components[i].pha == maximum) {
      amaxidx[nmaxidx++] = i;
    }
    if (components[i].pha == minimum) {
      aminidx[nminidx++] = i;
    }
  }

  // Case with 2 events have same energy, pick first in list.
  int equalenergyfound = 0;
  long maxidx=0, minidx=0, mididx=0;
  maxidx = min(amaxidx, nmaxidx);
  if (nmaxidx > 1) {
    mididx = max(amaxidx, nmaxidx);
    equalenergyfound = 1;
  }
  minidx = min(aminidx, nminidx);
  if (nminidx > 1) {
    mididx = max(aminidx, nminidx);
    equalenergyfound = 1;
  }

  // Check if event with maximum PHA is our main event:
  if ((event.xi == components[maxidx].xi) && (event.yi == components[maxidx].yi)) {
    
    // We found doubles:
    if (2 == event.patnum) {
      if ((1 == equalenergyfound) && (minidx == maxidx)) {
	minidx = 1;
	// Because we have two components and the event with the maximum energy
	// was already checked for to be the components[maxidx].
	// For equalenergyfound maxidx is the first in the amaxidx list,
	// therefore the minidx MUST be the other one (1).
      }

      // Check position of second event.
      if (event.yi > components[minidx].yi) event.patid = 1;
      if (event.xi < components[minidx].xi) event.patid = 2;
      if (event.yi < components[minidx].yi) event.patid = 3;
      if (event.xi > components[minidx].xi) event.patid = 4;

      // Case of bad corner event, xi and yi differ.
      if ((event.xi != components[minidx].xi) && (event.yi != components[minidx].yi))
	event.patid = -2;
    } // END if double.

    // Triples:
    if (3 == event.patnum) {
      // Get index of 3rd element.
      if (0 == equalenergyfound) {
	long amididx[ARRAY_LENGTH], nmididx=0;
	for (i=0; i<ncomponents; i++) {
	  if ((components[i].pha < maximum) && (components[i].pha > minimum)) {
	    amididx[nmididx++] = i;
	  }
	}
	// Consistency check:
	assert (1 == nmididx);
	mididx = amididx[0];
      } // Else use setting found above.

      // Search for position of element with minimum energy with respect to event
      // with maximum energy. Then find position  of third element.
      if ((event.yi-1 == components[minidx].yi) && (event.xi == components[minidx].xi)) {
	if ((event.xi-1 == components[mididx].xi) && (event.yi == components[mididx].yi))
	  event.patid = 5;
	if ((event.xi+1 == components[mididx].xi) && (event.yi == components[mididx].yi))
	  event.patid = 6;
      }

      if ((event.yi+1 == components[minidx].yi) && (event.xi == components[minidx].xi)) {
	if ((event.xi-1 == components[mididx].xi) && (event.yi == components[mididx].yi))
	  event.patid = 8;
	if ((event.xi+1 == components[mididx].xi) && (event.yi == components[mididx].yi))
	  event.patid = 7;
      }

      if ((event.xi-1 == components[minidx].xi) && (event.yi == components[minidx].yi)) {
	if ((event.yi-1 == components[mididx].yi) && (event.xi == components[mididx].xi))
	  event.patid = 5;
	if ((event.yi+1 == components[mididx].yi) && (event.xi == components[mididx].xi))
	  event.patid = 8;
      }

      if ((event.xi+1 == components[minidx].xi) && (event.yi == components[minidx].yi)) {
	if ((event.yi-1 == components[mididx].yi) && (event.xi == components[mididx].xi))
	  event.patid = 6;
	if ((event.yi+1 == components[mididx].yi) && (event.xi == components[mididx].xi))
	  event.patid = 7;
      }

      // If event build out of 3 events doesn't have patid up to now,
      // it is not a valid triple.
      if (-1 == event.patid) event.patid = 3000;
    } // END triples.

    // Quadruples:
    if (4 == event.patnum) {
      // Check if events build a 2x2 matrix.

      // A quad has four elements; here we determine the elements which are not
      // min or max elements of the components list. !UGLY!
      long extreme[2] = {minidx, maxidx};
      // Index of the central components:
      long centralidx[4], ncentralidx=0;
      for (i=0; i<4; i++) {
	if ((i != extreme[0]) && (i != extreme[1])) {
	  centralidx[ncentralidx++] = i;
	}
      }
      assert(2<=ncentralidx); 
      // TODO: quadruples with equal PHA values in all 4 pixels are
      // marked as invalid!!
      int pass = 0;

      // Now we test if the central elements share column (xi) and row (yi) 
      // with the min / max events. Both criteria must be satisfied.
      // The central events are on the corners of a 2x2 matrix.
      for (i=0; i<2; i++) {
	if (((components[centralidx[i]].yi == components[minidx].yi) &&
	     (components[centralidx[i]].xi == components[maxidx].xi)) ||
	    ((components[centralidx[i]].yi == components[maxidx].yi) &&
	     (components[centralidx[i]].xi == components[minidx].xi))) {
	  pass++;
	}
      }

      // If pass is 2 it has the valid shape of a quad.
      // Now determine the orientation of the pattern.
      if (2 == pass) {
	// Position of events with minimum and maximum energy must be over corner.
	// If this is true, set patid number.
	if ((event.yi-1 == components[minidx].yi) && 
	    (event.xi-1 == components[minidx].xi)) event.patid = 9;
	if ((event.yi-1 == components[minidx].yi) && 
	    (event.xi+1 == components[minidx].xi)) event.patid = 10;
	if ((event.yi+1 == components[minidx].yi) && 
	    (event.xi+1 == components[minidx].xi)) event.patid = 11;
	if ((event.yi+1 == components[minidx].yi) && 
	    (event.xi-1 == components[minidx].xi)) event.patid = 12;
      }

      // If event build out of 4 events doesn't have a patid up to now,
      // it is not a valid quadruple.
      if (-1 == event.patid) event.patid = 4000;
    } // END if quadruple

  } // END check if event with maximum energy is main event.
      
  return(event);
} // END pattern_id()



int wfi_pattern_recombination_main() {
  struct Parameters parameters;
  WFIEventFile eventfile;
  WFIEventFile patternfile;
  WFIEvent** ccdarr=NULL;

  char msg[MAXMSG];
  int status = EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("wfi_pattern_recombination");
  set_toolversion("0.01");

  do { // ERROR handling loop

    // Read parameters by PIL:
    status = pattern_recombination_getpar(&parameters);
    if (EXIT_SUCCESS!=status) break;

    // Initialize HEADAS random number generator.
    HDmtInit(1);

    // Read the EBOUNDS from the detector response file.
    struct RMF* rmf = loadRMF(parameters.response_filename, &status);
    if (EXIT_SUCCESS!=status) break;

    // Open the INPUT event file:
    status = openWFIEventFile(&eventfile, parameters.eventlist_filename, READONLY);
    if (EXIT_SUCCESS!=status) break;

    // Check if the input file is empty:
    if (0 == eventfile.rows) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: input file '%s' does not contain any events!\n", 
	      parameters.eventlist_filename);
      HD_ERROR_THROW(msg, status);
      break;
    }

    // Create a new OUTPUT event / pattern file:
    // Set the event list template file for the different WFI modes:
    char template_filename[MAXMSG];
    strcpy(template_filename, parameters.templatedir);
    strcat(template_filename, "/");
    if (16==eventfile.columns) {
      strcat(template_filename, "wfi.window16.eventlist.tpl");
    } else if (1024==eventfile.columns) {
      strcat(template_filename, "wfi.full1024.eventlist.tpl");
    } else if (64==eventfile.columns) {
      strcat(template_filename, "wfi.labor64.eventlist.tpl");
    } else {
      status = EXIT_FAILURE;
      sprintf(msg, "Error: detector width (%d pixels) is not supported!\n", eventfile.columns);
      HD_ERROR_THROW(msg, status);
      break;
    }
    // Create and open the new file
    status = openNewWFIEventFile(&patternfile, parameters.pattern_filename, template_filename);
    if (EXIT_SUCCESS!=status) break;


    // Get memory for the CCD array.
    ccdarr = (WFIEvent**)malloc(eventfile.columns*sizeof(WFIEvent*));
    if (NULL==ccdarr) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation failed!\n", status);
      break;
    } 
    long i;
    for (i=0; i<eventfile.columns; i++) {
      ccdarr[i] = (WFIEvent*)malloc(eventfile.rows*sizeof(WFIEvent));
      if (NULL==ccdarr[i]) {
	status = EXIT_FAILURE;
	HD_ERROR_THROW("Error: memory allocation failed!\n", status);
	break;
      }
    }
    if (EXIT_SUCCESS!=status) break;
    

    // Starting values for the emptyevent in the pattern search.
    WFIEvent emptyevent = {
      .pha = -1,
      .xi = -1,
      .yi = -1,
      .frame = -1,
      .patnum = -1,
      .patid = 0,
      .pileup = 0
    };


    // List of events belonging to the same frame.
    WFIEvent eventlist[ARRAY_LENGTH];
    // Components of a split pattern.
    WFIEvent components[ARRAY_LENGTH];
    // Read the first event from the input event file.
    long neventlist=1;
    status=WFIEventFile_getNextRow(&eventfile, &(eventlist[0]));
    if(EXIT_SUCCESS!=status) break;
    
    // Read in all subsequent events from the event file.
    WFIEvent newevent;
    while ((EXIT_SUCCESS==status) && (0==EventFileEOF(&eventfile.generic))) {

      // Read the next event from the FITS file.
      status=WFIEventFile_getNextRow(&eventfile, &newevent);
      if(EXIT_SUCCESS!=status) break;

      // Check if the new event belongs to the same frame as the previous ones:
      if (newevent.frame != eventlist[0].frame) {

	// Perform pattern recognition on this frame.
	// Write current events onto WFI detector pixel matrix:

	// Clear CCD array
	long j;
	for (i=0; i<eventfile.columns; i++) {
	  for (j=0; j<eventfile.rows; j++) {
	    ccdarr[i][j] = emptyevent;
	  }
	}

	assert(neventlist<ARRAY_LENGTH);
	for (i=0; i<neventlist; i++) {
	  ccdarr[eventlist[i].xi][eventlist[i].yi] = eventlist[i];
	  ccdarr[eventlist[i].xi][eventlist[i].yi].patnum = 0;
	}
	
	// Find singles and multiples:
	for (i=0; (i<neventlist)&&(EXIT_SUCCESS==status); i++) {
	  // Combine event at this position:
	  int column = eventlist[i].xi;
	  int row    = eventlist[i].yi;
	  // Check event, if it has not yet been dealt with:
	  if (-1!=ccdarr[column][row].patnum) { // Choose only pixels with events.
	    // Initialize maxevent (holding the recombined photon energy at the position
	    // of the event with the maximum photon energy.
	    WFIEvent maxevent = emptyevent;
	    // Becomes a list, containing the events the current pattern is build from.
	    long ncomponents = 0;

	    WFIEvent event = mark(row, column, eventfile.rows, eventfile.columns,
				  ccdarr, &maxevent, components, &ncomponents, rmf);

	    // Call pattern_id() to check for valild patterns and label them.
	    // Not needed for border events and singles.
	    // Here the generated pattern is labeled.
	    WFIEvent pattern;
	    if (1 == event.patnum) {
	      event.patid = 0; // Set single patid.
	      pattern = event;
	    } else {
	      pattern = pattern_id(components, ncomponents, event);
	    }

	    // Write the new pattern to the OUTPUT pattern file.
	    status = addWFIEvent2File(&patternfile, &pattern);
	    if (EXIT_SUCCESS!=status) break;

	  } // END choose only pixels with events.
	} // END of loop over all photons in this frame.

	//---------
	// Start new frame:
	eventlist[0] = newevent;
	neventlist = 1;

      } else {
	// Event belongs to the same frame as the previous ones, therefore just 
	// add it to the existing list.
	eventlist[neventlist++] = newevent;
      }
      
    } // END of loop over all events in the event file.
            
  } while(0); // End of error handling loop


  // --- Clean Up ---

  // Free memory from the CCD array.
  if (NULL!=ccdarr) {
    long i;
    for (i=0; i<eventfile.columns; i++) {
      if (NULL!=ccdarr[i]) {
	free(ccdarr[i]);
      }
    }
    free(ccdarr);
  }
      
  // Close the event files:
  closeWFIEventFile(&eventfile);
  closeWFIEventFile(&patternfile);

  // Release HEADAS random number generator.
  HDmtFree();

  return(status);
}


