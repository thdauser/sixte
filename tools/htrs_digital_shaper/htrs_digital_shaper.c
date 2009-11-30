#include "htrs_digital_shaper.h"


int htrs_digital_shaper_main() {
  struct Parameters parameters;

  // Input event file containing all events.
  HTRSEventFile input_eventfile;
  // Output event file containing only events that are properly shaped.
  HTRSEventFile output_eventfile;
  // Currently regarded event.
  HTRSEvent event;
  // Buffer for scanning the event file.
  HTRSEvent eventbuffer;

  // Period before and after each event that may not contain
  // additional events (in the same pixel) ([s]).
  double empty_period;

  int status = EXIT_SUCCESS;

  // Register HEATOOL
  set_toolname("htrs_digital_shaper");
  set_toolversion("0.01");

  do { // ERROR handling loop

    // Read parameters by PIL:
    status = htrs_digital_shaper_getpar(&parameters);
    if (EXIT_SUCCESS!=status) break;

    empty_period = 1./parameters.frequency * parameters.nsamplings;
    headas_chat(1, "HTRS digital shaper (sampling frequency %.2lf MHz) ...\n",
		parameters.frequency/1.e6);
    headas_chat(1, "minimum time between 2 subsequent events: %lf us\n", empty_period*1.e6);

    // Open the input event file.
    status=openHTRSEventFile(&input_eventfile, parameters.input_eventlist_filename, READWRITE);
    if (EXIT_SUCCESS!=status) break;
    // Create and open the output event file.
    status=openNewHTRSEventFile(&output_eventfile, parameters.output_eventlist_filename,
				parameters.eventlist_template);
    if (EXIT_SUCCESS!=status) break;


    // Loop over all events in the event file.
    while((EXIT_SUCCESS==status) && (0==EventFileEOF(&input_eventfile.generic))) {
      
      // Read the next event from the FITS file.
      status=HTRSEventFile_getNextRow(&input_eventfile, &event);
      if(EXIT_SUCCESS!=status) break;

      // Check the events before and after the current one within the specified 
      // time span.
      int properly_shaped=1;
      // Former events:
      long row = input_eventfile.generic.row-1;
      while ((1==EventFileRowIsValid(&input_eventfile.generic, row)) &&
	     (1==properly_shaped)){
	status = HTRSEventFile_getRow(&input_eventfile, &eventbuffer, row);
	if (EXIT_SUCCESS!=status) break;
	if (event.time - eventbuffer.time > empty_period) break;
	if (event.pixel == eventbuffer.pixel) properly_shaped=0;
	row--;
      }
      if (EXIT_SUCCESS!=status) break;
      // Subsequent events:
      row = input_eventfile.generic.row + 1;
      while ((1==EventFileRowIsValid(&input_eventfile.generic, row)) &&
	     (1==properly_shaped)){
	status = HTRSEventFile_getRow(&input_eventfile, &eventbuffer, row);
	if (EXIT_SUCCESS!=status) break;
	if (eventbuffer.time - event.time > empty_period) break;
	if (event.pixel == eventbuffer.pixel) properly_shaped=0;
	row++;
      }
      if (EXIT_SUCCESS!=status) break;

      // If the event was properly shaped, add it to the output event file.
      if (1==properly_shaped) {
	status=addHTRSEvent2File(&output_eventfile, &event);
	if (EXIT_SUCCESS!=status) break;
      }
      
    } // End of loop over all events in the event file
    
  } while(0); // End of error handling loop

  // --- Clean Up ---

  // Close the event files.
  closeHTRSEventFile(&input_eventfile);
  closeHTRSEventFile(&output_eventfile);

  return(status);
}



int htrs_digital_shaper_getpar(struct Parameters* parameters)
{
  int status = EXIT_SUCCESS;

  if ((status = PILGetFname("input_eventlist_filename", parameters->input_eventlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the input event file!\n", status);
  }

  else if ((status = PILGetFname("output_eventlist_filename", parameters->output_eventlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the output event file!\n", status);
  }

  else if ((status = PILGetReal("frequency", &parameters->frequency))) {
    HD_ERROR_THROW("Error reading the sampling frequency!\n", status);
  }

  else if ((status = PILGetInt("nsamplings", &parameters->nsamplings))) {
    HD_ERROR_THROW("Error reading the number of required samplings!\n", status);
  }

  // Get the name of the FITS template directory.
  // First try to read it from the environment variable.
  // If the variable does not exist, read it from the PIL.
  else { 
    char* buffer;
    if (NULL!=(buffer=getenv("SIXT_FITS_TEMPLATES"))) {
      strcpy(parameters->eventlist_template, buffer);
    } else {
      if ((status = PILGetFname("fits_templates", parameters->eventlist_template))) {
	HD_ERROR_THROW("Error reading the path of the FITS templates!\n", status);
      }
    }
  }
  if (EXIT_SUCCESS!=status) return(status);
  // Set the event list template file:
  strcat(parameters->eventlist_template, "/htrs.eventlist.tpl");

  return(status);
}





