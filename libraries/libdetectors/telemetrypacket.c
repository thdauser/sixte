#include "telemetrypacket.h"


int initTelemetryPacket(TelemetryPacket* packet, int nbits)
{
  // Allocate memory for the TelemetryPacket data.
  packet->data = (unsigned char*)malloc(nbits*sizeof(unsigned char));
  if (NULL==packet->data) {
    HD_ERROR_THROW("Error: could not allocate memory for TelemetryPacket!\n", 
		   EXIT_FAILURE);
    return(EXIT_FAILURE);
  }

  // Set the total number of bits and the pointer to the current bit.
  packet->nbits = nbits;
  packet->current_bit=0;  

  return(EXIT_SUCCESS);
}



void cleanupTelemetryPacket(TelemetryPacket* packet)
{
  if (packet->nbits>0) {
    if (NULL!=packet->data) {
      free(packet->data);
    }
  }
  packet->nbits=0;
  packet->current_bit=0;
}



void newTelemetryPacket(TelemetryPacket* packet)
{
  int count;
  for(count=0; count<packet->nbits; count++) {
    packet->data[count]=0;
  }
  packet->current_bit=0;
}


