#include "htrstelemetrypacket.h"


int initHTRSTelemetryPacket(HTRSTelemetryPacket* packet,
			    struct HTRSTelemetryPacketParameters* parameters)
{
  packet->n_header_bits = parameters->n_header_bits;

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



int addData2HTRSTelemetryPacket(HTRSTelemetryPacket* packet, unsigned char* data, int nbits)
{
  return(addData2TelemetryPacket(&packet->tp, data, nbits));
}



int writeHTRSTelemetryPacket2File(HTRSTelemetryPacket* packet, FILE* file)
{
  return(writeTelemetryPacket2File(&packet->tp, file));
}


