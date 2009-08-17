#include "analyse_xms_events.h"


int analyse_xms_events_main() {
  struct Parameters parameters;
  XMSEventFile eventfile;
  int status = EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("scan_tes_events");
  set_toolversion("0.01");


  do { // ERROR handling loop

    // Read parameters by PIL:
    status = analyse_xms_events_getpar(&parameters);
    if (EXIT_SUCCESS!=status) break;

    // Open the event file.
    status=openXMSEventFile(&eventfile, parameters.eventlist_filename, READWRITE);
    if (EXIT_SUCCESS!=status) break;

  } while(0); // End of error handling loop


  // --- clean up ---
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



void insert_tes_event(double* event_list, double time)
{
  int t;

  // Shift the former events in the list of this pixel:
  for (t=TES_EVENT_LIST_ENTRIES-1; t>0; t--) {
    event_list[t] = event_list[t-1];
  }
  // Store the new event in the list:
  event_list[0] = time;
}



int analyse_xms_events_main() {
  struct Event_List_File event_list_file;
  double*** detector=NULL;
  int det_width;

  double t_start, t_end; // start and end time of observation (FITS header keywords)
  double time_unit;      // detector intrinsic time unit
  int units_before_pulse, units_after_pulse;
  double time_before_pulse, time_after_pulse;
  double equilibration_time;

  // The observation is performed only in an interval from 
  // ( t_start + equilibration_time  ;  t_end - equilibration_time )
  // in order to avoid border effects.
  // This interval is stored in these variables.
  double t_start_obs, t_end_obs;

  long event_type[3] = {0, 0, 0};

  char msg[MAXMSG];             // error output buffer
  int status = EXIT_SUCCESS;


  // Register HEATOOL
  set_toolname("scan_tes_events");
  set_toolversion("0.01");


  do { // ERROR handling loop
    // Read parameters by PIL:
    status = scan_tes_events_getpar(event_list_file.filename, &time_unit,
				    &units_before_pulse, &units_after_pulse);

    // If an error has occurred during parameter input, 
    // break the program execution.
    if (status!=EXIT_SUCCESS) break;
      
    // Determine parameters from input data:
    time_before_pulse = units_before_pulse * time_unit;
    time_after_pulse  = units_after_pulse  * time_unit;
    equilibration_time = time_before_pulse + time_after_pulse;
    headas_chat(5, "Time before high precision event: %lf s\n", time_before_pulse);
    headas_chat(5, "Time after high precision event: %lf s\n", time_after_pulse);
    headas_chat(5, "Equilibration time: %lf s\n", equilibration_time);

  
    // Open the input FITS file to get the event list table:
    open_event_list_file(&event_list_file, &status);

    // Get the width of the detector from the header keywords.
    char comment[MAXMSG];   // input buffer for header comments
    if (fits_read_key(event_list_file.fptr, TINT, "DETWIDTH", &det_width, 
		      comment, &status)) break;
    if (fits_read_key(event_list_file.fptr, TDOUBLE, "TSTART", &t_start, 
		      comment, &status)) break;
    if (fits_read_key(event_list_file.fptr, TDOUBLE, "TEND", &t_end, 
		      comment, &status)) break;
    t_start_obs = t_start + equilibration_time;
    t_end_obs   = t_end   - equilibration_time;
    headas_chat(5, "Observation time interval: from %lf s to %lf s\n", 
		t_start_obs, t_end_obs);


    // Get memory for the event lists in each individual pixel of the detector.
    int x, y, t;
    detector = (double ***) malloc (det_width * sizeof(double **));
    if (detector) {
      for (x = 0; x < det_width; x++) {
	detector[x] = (double **) malloc(det_width * sizeof(double *));
	if (detector[x]) {
	  for (y = 0; y < det_width; y++) {
	    detector[x][y] = 
	      (double *) malloc(TES_EVENT_LIST_ENTRIES * sizeof(double));
	    if (detector[x][y]) {
	      // clear the detector array
	      for (t=0; t<TES_EVENT_LIST_ENTRIES; t++) {
		detector[x][y][t] = 0.;
	      }
	    } else {
	      status=EXIT_FAILURE;
	      sprintf(msg, "Error: not enough memory!\n");
	      HD_ERROR_THROW(msg,status);
	      break;
	    }
	  }
	} else {
	  status=EXIT_FAILURE;
	  sprintf(msg, "Error: not enough memory!\n");
	  HD_ERROR_THROW(msg,status);
	  break;
	}
      }
    } else {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }  // END of memory allocation for detector
    


    // Read the event list and perform TES event analysis
    // (high precision or degraded event quality etc.).
    headas_chat(5, "processing events ...\n");
    struct Event event;  // input buffer
    for(event_list_file.row=0; 
	(event_list_file.row<event_list_file.nrows)&&(status==EXIT_SUCCESS); 
	event_list_file.row++) {
      
      if(get_eventtbl_row(event_list_file, &event, &status)) break;


      // Insert the new event in the detector pixel list.
      insert_tes_event(detector[event.xi][event.yi], event.time);


      // Analyse the event list of the affected pixel.
      analyse_tes_event(detector[event.xi][event.yi],
			t_start+equilibration_time, t_end-equilibration_time,
			time_before_pulse, time_after_pulse, event_type);

    }  // END of scanning the event list


    // Work on the events that remain in the detector pixel event list.
    for (x=0; x<det_width; x++) {
      for (y=0; y<det_width; y++) {
	for (t=0; t<2; t++) {
	  // Insert dummy photon at the end of the event list time interval.
	  insert_tes_event(detector[x][y], t_end);
	
	  // Analyse pixel.
	  analyse_tes_event(detector[x][y],
			    t_start+equilibration_time, t_end-equilibration_time,
			    time_before_pulse, time_after_pulse, event_type);
	}
      }
    }

	


    // Display the results on STDOUT:
    headas_chat(0, "Number of total events: %ld\n", event_type[0]);
    headas_chat(0, "Number of high precision events: %ld\n", event_type[1]);
    headas_chat(0, "Number of high AND degraded precision events: %ld\n", 
		event_type[2]);
    headas_chat(0, "Observed time interval: %lf s\n", 
		t_end - t_start - 2*equilibration_time);

  } while(0);  // END of error handling loop

  if (event_list_file.fptr != NULL) fits_close_file(event_list_file.fptr, &status);

  return(status);
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





