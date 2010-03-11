#include "htrs_digital_shaper.h"


int htrs_digital_shaper_main() {
  struct Parameters parameters;

  // Input event file containing all events.
  HTRSEventFile input_eventfile;
  // Output event file containing only events that are properly shaped.
  HTRSEventFile output_eventfile;

  // Currently processed event.
  HTRSEvent event;
  // Buffer for scanning the event file.
  HTRSEvent eventbuffer;

  // Period before and after each event that may not contain
  // additional events (in the same pixel) ([s]).
  double shaping_time;

  // ADC history.
  // For each ADC the last NTIMEBINS channel values are stored in
  // order to apply the different filters.
  long ADC[NCHANNELS];
  // End of the interval of the next reset.
  double reset_time[NCHANNELS];
  // Channel counter.
  int channel;

  // Width of an ADC channel in [keV].
  double adc_channel_width;

  int status = EXIT_SUCCESS;

  // Register HEATOOL
  set_toolname("htrs_digital_shaper");
  set_toolversion("0.01");

  do { // ERROR handling loop

    // Read parameters by PIL:
    status = htrs_digital_shaper_getpar(&parameters);
    if (EXIT_SUCCESS!=status) break;

    // Calculate the time step between two subsequent ADC calls.
    headas_chat(1, "HTRS digital shaper (sampling frequency %.2lf MHz, "
		"%d samplings) ...\n",
		parameters.frequency/1.e6,
		parameters.nsamplings);
    shaping_time = 1./parameters.frequency * parameters.nsamplings;
    headas_chat(1, " -> shaping time: %lf mu s\n", shaping_time*1.e6);

    // Calculate the ADC properties according to the current settings.
    // One ADC channel corresponds to ??? keV:
    adc_channel_width = adc_input_range / 
      (sdd_output * parameters.preamp_gain * adc_N_channels);

    // Initialize arrays.
    for (channel=0; channel<NCHANNELS; channel++) {
      ADC[channel] = 0;
      reset_time[channel] = 0.;
    }

    // Open the input event file.
    status=openHTRSEventFile(&input_eventfile, 
			     parameters.input_eventlist_filename, 
			     READWRITE);
    if (EXIT_SUCCESS!=status) break;

    // Create and open the output event file.
    status=openNewHTRSEventFile(&output_eventfile, 
				parameters.output_eventlist_filename,
				parameters.eventlist_template);
    if (EXIT_SUCCESS!=status) break;

    
    // Loop over all events in the event file.
    while((EXIT_SUCCESS==status) && (0==EventFileEOF(&input_eventfile.generic))) {
      
      // Read the next event from the FITS file.
      status=HTRSEventFile_getNextRow(&input_eventfile, &event);
      if(EXIT_SUCCESS!=status) break;

      int properly_shaped=1;

      // Check if the pixel is current in a reset interval. In that
      // case the event must be discarded.
      if (event.time < reset_time[event.pixel]) properly_shaped=0;

      // Check the events before and after the current one within the specified 
      // time span.
      // Former events:
      long row = input_eventfile.generic.row-1;
      while ((1==EventFileRowIsValid(&input_eventfile.generic, row)) &&
	     (1==properly_shaped)){
	status = HTRSEventFile_getRow(&input_eventfile, &eventbuffer, row);
	if (EXIT_SUCCESS!=status) break;
	if (event.time - eventbuffer.time > shaping_time) break;
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
	if (eventbuffer.time - event.time > shaping_time) break;
	if (event.pixel == eventbuffer.pixel) properly_shaped=0;
	row++;
      }
      if (EXIT_SUCCESS!=status) break;

      // If the event was properly shaped, add it to the output event file.
      // TODO Add the event to the output list in any case, 
      // but specify an addition column with the type of the event 
      // (properly detected, during reset, during dead time, ...).
      if (1==properly_shaped) {
	status=addHTRSEvent2File(&output_eventfile, &event);
	if (EXIT_SUCCESS!=status) break;
      }

      // Increase the ADC channel by the value corresponding to 
      // the photon energy.
      ADC[event.pixel] += event.energy / adc_channel_width;
      // (both quantities are in [keV].)
      
      // Check if a reset must be applied to the SDD.
      if (ADC[event.pixel] >= adc_reset_threshold) {
	ADC[event.pixel] = 0;
	reset_time[event.pixel] = event.time + reset_duration;
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

  else if ((status = PILGetFname("output_eventlist_filename", 
				 parameters->output_eventlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the output event file!\n", status);
  }

  else if ((status = PILGetReal("frequency", &parameters->frequency))) {
    HD_ERROR_THROW("Error reading the sampling frequency!\n", status);
  }

  else if ((status = PILGetReal("preamp_gain", &parameters->preamp_gain))) {
    HD_ERROR_THROW("Error reading the gain of the preamplifier!\n", status);
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





