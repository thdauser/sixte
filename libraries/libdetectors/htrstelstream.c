#include "htrstelstream.h"


HTRSTelStream* getHTRSTelStream(struct HTRSTelStreamParameters* parameters, 
				int* status)
{
  HTRSTelStream* htstream = (HTRSTelStream*)malloc(sizeof(HTRSTelStream));
  if (NULL==htstream) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for HTRSTelStream!\n",
		   *status);
    return(htstream);
  }
  
  // Set the object properties.
  htstream->n_header_bits = parameters->n_header_bits;
  htstream->n_bin_bits    = parameters->n_bin_bits;
  htstream->n_channels    = parameters->n_channels;
  htstream->n_bins        = parameters->n_bins;
  htstream->spectrum_time = 0.;
  htstream->integration_time = parameters->integration_time;
  htstream->tp            = NULL;
  htstream->output_file   = NULL;
  htstream->spectrum      = NULL;
  htstream->chans2bins    = parameters->chans2bins;

  // Determine the maximum number of counts per spectrum bin.
  // This number is determined by the number of available bits per bin.
  assert(htstream->n_bin_bits <= 8);
  htstream->max_counts = 0x0001 << (htstream->n_bin_bits-1);

  // Get a new TelemetryPacket object.
  htstream->tp = getTelemetryPacket(parameters->n_packet_bits, status);
  if(EXIT_SUCCESS!=*status) return(htstream);
  // Start a new TelemetryPacket.
  newHTRSTelStreamPacket(htstream);

  // Allocate memory for the binned spectrum.
  htstream->spectrum=(int*)malloc(htstream->n_bin_bits*sizeof(int));
  if (NULL==htstream->spectrum) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for HTRSTelStream!\n", 
		   *status);
    return(htstream);
  }

  // Open the binary output file.
  htstream->output_file = fopen(parameters->output_filename, "w+");
  if (NULL==htstream->output_file) {
    char msg[MAXMSG];
    *status=EXIT_FAILURE;
    sprintf(msg, "Error: output file '%s' could not be opened!\n", 
	    parameters->output_filename);
    HD_ERROR_THROW(msg, *status);
    return(htstream);
  }
  
  return(htstream);
}



void freeHTRSTelStream(HTRSTelStream* stream)
{
  if(NULL!=stream) {
    // Call the destructor of the TelemetryPacket object.
    freeTelemetryPacket(stream->tp);
    // Free the binned spectrum.
    if (NULL!=stream->spectrum) free(stream->spectrum);
    // Close the binary output file.
    fclose(stream->output_file);
    // Free the memory of the HTRSTelStream object itself.
    free(stream);
  }
}



int addEvent2HTRSTelStream(HTRSTelStream* stream, HTRSEvent* event)
{
  int status=EXIT_SUCCESS;

  // Check if the event can be added to the existing spectrum, or whether
  // a new telemetry packet has to be started. In the latter case the old
  // packet will be written to the binary output file.
  while (event->time > stream->spectrum_time+stream->integration_time) {
    // A new spectrum has to be started.

    // Check whether the TelemetryPacket is full. In that case write the
    // content to the binary output file and start a new TelemetryPacket.
    if (availableBitsInTelemetryPacket(stream->tp) < 
	stream->n_bin_bits*stream->n_bins) {

      // Write the content of the packet to the binary output file.
      status= writeTelemetryPacket2File(stream->tp, stream->output_file);
      if (EXIT_SUCCESS!=status) break;

      // Start a new packet.
      newHTRSTelStreamPacket(stream);
    }

    // TODO Add the binned spectrum to the TelemetryPacket.

    // Clear the spectrum.
    int count;
    for (count=0; count<stream->n_bins; count++) {
      stream->spectrum[count] = 0;
    }

    // Update the spectrum time.
    stream->spectrum_time += stream->integration_time;
  }
  if (EXIT_SUCCESS!=status) return(status);

  // The event can be added to the current spectrum.
  assert(event->pha>0); assert(event->pha<=stream->n_channels);
  stream->spectrum[stream->chans2bins[event->pha-1]]++;

  return(status);
}



int completeHTRSTelStream(HTRSTelStream* stream)
{
  // TODO

  return(EXIT_SUCCESS);
}



void newHTRSTelStreamPacket(HTRSTelStream* stream)
{
  // Initialize a new TelemetryPacket.
  newTelemetryPacket(stream->tp);

  // TODO Insert the header a the beginning of the new packet.

}


