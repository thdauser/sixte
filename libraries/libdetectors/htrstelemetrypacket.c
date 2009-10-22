#include "htrstelemetrypacket.h"


int initHTRSTelemetryPacket(HTRSTelemetryPacket* packet,
			    struct HTRSTelemetryPacketParameters* parameters)
{
  // The total number of bits per packet must be a multiple of 8.
  assert(parameters->n_packet_bits % 8 == 0);

  packet->n_header_bits   = parameters->n_header_bits;

  // Currently the total number of bits per spectrum must be a multiple of 8, 
  // i.e., we are only dealing with complete bytes.
  packet->n_spectrum_bits = parameters->n_spectrum_bits;
  assert(packet->n_spectrum_bits % 8 == 0);

  packet->n_bin_bits      = parameters->n_bin_bits;

  // Determine the maximum number of counts per spectral bin.
  unsigned char maximum=0;
  int count;
  for (count=0; count<packet->n_bin_bits; count++) {
    maximum = (maximum<<1) + 1;
  }
  packet->max_counts = (int)maximum;

  // Allocate memory for the byte buffer.
  packet->byte_buffer=(unsigned char*)malloc(parameters->n_packet_bits/8*sizeof(unsigned char*));
  if (NULL==packet->byte_buffer) {
    HD_ERROR_THROW("Error: Could not allocate memory for byte buffer!\n", EXIT_FAILURE);
    return(EXIT_FAILURE);
  }

  return(initTelemetryPacket(&packet->tp, parameters->n_packet_bits));
}



void cleanupHTRSTelemetryPacket(HTRSTelemetryPacket* packet)
{
  cleanupTelemetryPacket(&packet->tp);
  if (NULL!=packet->byte_buffer) free(packet->byte_buffer);
}



int newHTRSTelemetryPacket(HTRSTelemetryPacket* packet)
{
  // Initialize new empty Packet.
  newTelemetryPacket(&packet->tp);

  // Add the Packet header.
  int count;
  for (count=0; count*8<packet->n_header_bits; count++) {
    packet->byte_buffer[count] = 0;
  }

  return(addData2TelemetryPacket(&packet->tp, packet->byte_buffer, count*8));
}



int availableBitsInHTRSTelemetryPacket(HTRSTelemetryPacket* packet)
{
  return(availableBitsInTelemetryPacket(&packet->tp));
}



int addSpectrum2HTRSTelemetryPacket(HTRSTelemetryPacket* packet, int* spectrum, 
				    int n_spectrum_bins)
{
  // Convert the spectrum to binary format and add it to the TelemetryPacket.
  int count, byte_index, bit_in_byte;
  for(byte_index=0; byte_index<packet->n_spectrum_bits/8; byte_index++) {
    packet->byte_buffer[byte_index]=0;
  }
  for(count=0; count<n_spectrum_bins; count++) {
    // Check whether the maximum number of counts in this spectral bin is exceeded.
    // (Overflows due to the limited number of bits per bin have to be avoided.)
    if(spectrum[count]>packet->max_counts) {
      printf("Warning: overflow (maximum number of counts per bin exceeded)!\n");
      spectrum[count]=packet->max_counts;
    }

    byte_index = (count*packet->n_bin_bits)/8;
    bit_in_byte= (count*packet->n_bin_bits)%8;
    // Check whether the current spectral bin fits within the current
    // byte, ...
    if (bit_in_byte <= 8-packet->n_bin_bits) {
      packet->byte_buffer[byte_index] += 
	(unsigned char)(spectrum[count]<<(8-packet->n_bin_bits-bit_in_byte));
    } else {
      // ... or whether it overlaps with 2 bytes.
      packet->byte_buffer[byte_index] += 
	(unsigned char)(spectrum[count]>>(packet->n_bin_bits+bit_in_byte-8));
      packet->byte_buffer[byte_index+1] += 
	(unsigned char)(spectrum[count]<<(16-packet->n_bin_bits-bit_in_byte));
    }

  }
 
  return(addData2TelemetryPacket(&packet->tp, packet->byte_buffer, packet->n_spectrum_bits));
}



int writeHTRSTelemetryPacket2File(HTRSTelemetryPacket* packet, FILE* file)
{
  return(writeTelemetryPacket2File(&packet->tp, file));
}


