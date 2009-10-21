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



int availableBitsInTelemetryPacket(TelemetryPacket* packet)
{
  return(packet->nbits-packet->current_bit);
}



int addData2TelemetryPacket(TelemetryPacket* packet, unsigned char* data, int nbits)
{
  // Currently only multiple numbers of 8 bits (1 byte) are allowed for the number
  // of bytes to be transferred.
  assert(nbits%8==0);

  int index;
  for(index=0; index*8<nbits; index++) {
    // Copy the data to the internal storage of the TelemetryPacket.
    packet->data[packet->current_bit/8] = data[index];
    // Increase the internal bit counter.
    packet->current_bit += 8;
  }

  return(EXIT_SUCCESS);
}



int writeTelemetryPacket2File(TelemetryPacket* packet, FILE* file)
{
  // Write the bytes from the internal storage to the file.
  int nbytes = fwrite (packet->data, 1, packet->nbits/8, file);
  if (nbytes < packet->nbits) {
    return(EXIT_FAILURE);
  } else {
    return(EXIT_SUCCESS);
  }
}


