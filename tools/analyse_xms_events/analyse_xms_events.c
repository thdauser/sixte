#include "analyse_xms_events.h"


int analyse_xms_events_main() {
  struct Parameters parameters;
  XMSEventFile eventfile;

  int status = EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("analyse_xms_events");
  set_toolversion("0.01");

  do { // ERROR handling loop

    // Read parameters by PIL:
    status = analyse_xms_events_getpar(&parameters);
    if (EXIT_SUCCESS!=status) break;

    // Open the event file.
    status=openXMSEventFile(&eventfile, parameters.eventlist_filename, READWRITE);
    if (EXIT_SUCCESS!=status) break;


    // Loop over all events in the event file.
    while((EXIT_SUCCESS==status) && (0==EventFileEOF(&eventfile.generic))) {
      
      // Read the next event from the FITS file.
      XMSEvent event;
      status=XMSEventFile_getNextRow(&eventfile, &event);
      if(EXIT_SUCCESS!=status) break;


      // Check the events before and after the current one within the specified 
      // time spans.
      XMSEvent eventbuffer;
      int nbefore=0, nafter=0;
      // Former events:
      long row = eventfile.generic.row-1;
      while (1==EventFileRowIsValid(&eventfile.generic, row)) {
	status = XMSEventFile_getRow(&eventfile, &eventbuffer, row);
	if (EXIT_SUCCESS!=status) break;
	if (event.time - eventbuffer.time > 
	    parameters.units_before_pulse * parameters.time_unit) break;
	if ((event.xi == eventbuffer.xi) && (event.yi == eventbuffer.yi)) nbefore++;
	if (nbefore > 2) break; // Avoid too many unnecessary loop runs.
	row--;
      }
      if (EXIT_SUCCESS!=status) break;
      // Subsequent events:
      row = eventfile.generic.row + 1;
      while (1==EventFileRowIsValid(&eventfile.generic, row)) {
	status = XMSEventFile_getRow(&eventfile, &eventbuffer, row);
	if (EXIT_SUCCESS!=status) break;
	if (eventbuffer.time - event.time > 
	    parameters.units_after_pulse * parameters.time_unit) break;
	if ((event.xi == eventbuffer.xi) && (event.yi == eventbuffer.yi)) nafter++;
	if (nafter > 2) break; // Avoid too many unnecessary loop runs.
	row++;
      }
      if (EXIT_SUCCESS!=status) break;


      // Determine the grade of the event.
      if (0 == nbefore + nafter) {
	event.grade = 1; // High precision
      } else if (1 == nbefore + nafter) {
	event.grade = 2; // Intermediate precision
      } else {
	event.grade = 1000; // Worse
      }


      // Write the event information to the event file.
      status = XMSEventFile_writeRow(&eventfile, &event, eventfile.generic.row);
      if (EXIT_SUCCESS!=status) break;

    } // End of loop over all events in the event file
    
  } while(0); // End of error handling loop


  // --- Clean Up ---

  // Close the event file.
  closeXMSEventFile(&eventfile);

  return(status);
}



/*
void analyse_tes_event(
		       double* event_list,
		       double t_0, double t_1,
		       double time_before, double time_after,
		       long* event_type
		       )
{
  // Check the event in the middle of the event list (index [2]).
  if ((event_list[2] > t_0) && (event_list[2] < t_1)) {
    // The event lies within the valid observation interval!

    // Count total events within the observation interval
    event_type[0]++;


    // Check for HIGH PRECISION events.
    if ((event_list[2]-event_list[3] >= time_before) &&
	// distance to former event
	(event_list[1]-event_list[2] >= time_after))
        // enough relaxation time
      { // valid high precision event!
	event_type[1]++;
      }

    
    // Check for HIGH precision AND DEGRADED precision events:
    if ((event_list[2]-event_list[4]>=time_before) &&
	(event_list[0]-event_list[2]>=time_after) &&
	(event_list[1]-event_list[3]>=time_before+time_after))
      { // valid event (either high or degraded precision) !
	event_type[2]++;
      }
  } // END of check for observation time interval
}
*/


int analyse_xms_events_getpar(struct Parameters* parameters)
{
  int status = EXIT_SUCCESS;

  if ((status = PILGetFname("eventlist_filename", parameters->eventlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the input file!\n", status);
  }

  else if ((status = PILGetReal("time_unit", &parameters->time_unit))) {
    HD_ERROR_THROW("Error reading the detector intrinsic time unit!\n", status);
  }

  else if ((status = PILGetInt("units_before_pulse", &parameters->units_before_pulse))) {
    HD_ERROR_THROW("Error reading the time units required before a "
		   "high precision event!\n", status);
  }

  else if ((status = PILGetInt("units_after_pulse", &parameters->units_after_pulse))) {
    HD_ERROR_THROW("Error reading the time units required after a "
		   "high precision event!\n", status);
  }

  return(status);
}





