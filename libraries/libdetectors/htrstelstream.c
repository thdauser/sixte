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
  htstream->max_counts = 0;
  int count;
  for (count=0; count<htstream->n_bin_bits; count++) {
    htstream->max_counts = (htstream->max_counts<<1) +1;
  }

  // Get a new TelemetryPacket object.
  htstream->tp = getTelemetryPacket(parameters->n_packet_bits, status);
  if(EXIT_SUCCESS!=*status) return(htstream);
  // Start a new TelemetryPacket.
  *status=newHTRSTelStreamPacket(htstream);
  if (EXIT_SUCCESS!=*status) return(htstream);

  // Allocate memory for the binned spectrum.
  htstream->spectrum=(int*)malloc(htstream->n_bins*sizeof(int));
  if (NULL==htstream->spectrum) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for HTRSTelStream!\n", 
		   *status);
    return(htstream);
  }
  // Clear the spectrum.
  for (count=0; count<htstream->n_bins; count++) {
    htstream->spectrum[count] = 0;
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

  do { // Beginning of ERROR Handling loop.

    // Check if the event can be added to the existing spectrum, or whether
    // a new telemetry packet has to be started. In the latter case the old
    // packet will be written to the binary output file.
    while (event->time > stream->spectrum_time+stream->integration_time) {
      // A new spectrum has to be started. So add the current spectrum
      // to the telemetry packet.
      status=HTRSTelStreamAddSpec2Packet(stream);
      if (EXIT_SUCCESS!=status) break;
    } // END of while loop.
    if (EXIT_SUCCESS!=status) break;

    // The event can be added to the current spectrum.
    assert(event->pha>0); assert(event->pha<=stream->n_channels);
    // Check for buffer overflows.
    if (stream->spectrum[stream->chans2bins[event->pha-1]]<stream->max_counts) {
      stream->spectrum[stream->chans2bins[event->pha-1]]++;
    } else {
      printf("### Warning: Buffer Overflow!\n");
      printf("    Maximum number of counts exceeded in bin %d (%d/%d).\n",
	     stream->chans2bins[event->pha-1],
	     stream->spectrum[stream->chans2bins[event->pha-1]]+1,
	     stream->max_counts);
    }

  } while(0); // END of Error Handling loop.

  return(status);
}



int finalizeHTRSTelStream(HTRSTelStream* stream)
{
  // Add the current spectrum to the telemetry packet.
  int status=HTRSTelStreamAddSpec2Packet(stream);
  if (EXIT_SUCCESS!=status) return(status);

  // If the current TelemetryPacket contains some data 
  // (except for the header), write it to the output file.
  if (stream->tp->current_bit > stream->n_header_bits) {
    // Write the content of the packet to the binary output file.
    status = writeTelemetryPacket2File(stream->tp, stream->output_file);
    if (EXIT_SUCCESS!=status) return(status);
  }
  
  return(status);
}



int HTRSTelStreamAddSpec2Packet(HTRSTelStream* stream)
{
  unsigned char* byte_buffer=NULL;
  int status=EXIT_SUCCESS;

  do { // Beginning of ERROR Handling loop.

    // Check whether the TelemetryPacket is full. In that case write the
    // content to the binary output file and start a new TelemetryPacket.
    if (availableBitsInTelemetryPacket(stream->tp) < 
	stream->n_bin_bits*stream->n_bins) {
      
      // Write the content of the packet to the binary output file.
      status = writeTelemetryPacket2File(stream->tp, stream->output_file);
      if (EXIT_SUCCESS!=status) break;
	
      // Start a new packet.
      status=newHTRSTelStreamPacket(stream);
      if (EXIT_SUCCESS!=status) break;
    }
      
    // Add the binned spectrum to the TelemetryPacket.
    int n_bytes = (stream->n_bins*stream->n_bin_bits)/8;
    if (0!=(stream->n_bins*stream->n_bin_bits)%8) {
      n_bytes+=1;
    }
    if (NULL==byte_buffer) {
      byte_buffer=(unsigned char*)malloc(n_bytes*sizeof(unsigned char));
    }
    if (NULL==byte_buffer) {
      status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: Memory allocation failed for byte buffer in HTRSTelStream!\n",
		     status);
      break;
    }
    // Clean the byte buffer.
    int count;
    for (count=0; count<n_bytes; count++) {
      byte_buffer[count] = 0;
    }
    // Fill the byte buffer with data from the spectrum.
    int index, modulus;
    unsigned char byte;
    for (count=0; count<stream->n_bins; count++) {
      byte = (unsigned char)stream->spectrum[count];
      index = (count*stream->n_bin_bits)/8;
      modulus = (count*stream->n_bin_bits)%8;
      if (modulus+stream->n_bin_bits<=8) {
	byte_buffer[index]   = byte_buffer[index]   | 
	  (byte<<(8-stream->n_bin_bits-modulus));
      } else {
	byte_buffer[index]   = byte_buffer[index]   | 
	  (byte>>(stream->n_bin_bits+modulus-8));
	byte_buffer[index+1] = byte_buffer[index+1] | 
	  (byte<<(16-stream->n_bin_bits-modulus));
      }
    }
    // Insert the byte buffer to the TelemtryPacket.
    status = addData2TelemetryPacket(stream->tp, byte_buffer, 
				       stream->n_bins*stream->n_bin_bits);
    if (EXIT_SUCCESS!=status) break;
    
    // Clear the spectrum.
    for (count=0; count<stream->n_bins; count++) {
      stream->spectrum[count] = 0;
    }
      
    // Update the spectrum time.
    stream->spectrum_time += stream->integration_time;

  } while(0); // END of Error Handling loop.

  // --- Clean up ---
  if (NULL!=byte_buffer) free(byte_buffer);

  return(status);
}



int newHTRSTelStreamPacket(HTRSTelStream* stream)
{
  // Initialize a new TelemetryPacket.
  newTelemetryPacket(stream->tp);

  // Insert the header a the beginning of the new packet.
  assert(140==stream->n_header_bits);
  unsigned char bytes[18];
  bytes[0]=0x01;
  bytes[1]=0x23;
  bytes[2]=0x45;
  bytes[3]=0x67;
  bytes[4]=0x89;
  bytes[5]=0xAB;
  bytes[6]=0xCD;
  bytes[7]=0xEF;
  bytes[8]=0x01;
  bytes[9]=0x23;
  bytes[10]=0x45;
  bytes[11]=0x67;
  bytes[12]=0x89;
  bytes[13]=0xAB;
  bytes[14]=0xCD;
  bytes[15]=0xEF;
  bytes[16]=0x01;
  bytes[17]=0x23;
  // Insert the byte buffer to the TelemtryPacket.
  return(addData2TelemetryPacket(stream->tp, bytes, stream->n_header_bits));
}


