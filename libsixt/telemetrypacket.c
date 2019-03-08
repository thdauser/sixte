/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "telemetrypacket.h"


TelemetryPacket* getTelemetryPacket(int nbits, int* status)
{
  assert(nbits%8 == 0);

  TelemetryPacket* packet = (TelemetryPacket*)malloc(sizeof(TelemetryPacket));
  if (NULL==packet) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: could not allocate memory for TelemetryPacket!\n",
		   *status);
    return(packet);
  }

  // Allocate memory for the TelemetryPacket data.
  packet->data = (unsigned char*)malloc(nbits*sizeof(unsigned char));
  if (NULL==packet->data) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: could not allocate memory for TelemetryPacket!\n",
		   *status);
    return(packet);
  }

  // Set the total number of bits and the pointer to the current bit.
  packet->nbits = nbits;
  packet->current_bit=0;

  return(packet);
}



void freeTelemetryPacket(TelemetryPacket* packet)
{
  if (NULL!=packet) {
    if (packet->nbits>0) {
      if (NULL!=packet->data) {
	free(packet->data);
      }
    }
    packet->nbits=0;
    packet->current_bit=0;
    free(packet);
  }
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



int addData2TelemetryPacket(TelemetryPacket* packet, unsigned char* data,
			    int nbits)
{
  unsigned char byte=0; // Buffer for the data.

  int index; // Current byte.
  unsigned char modulus; // Current bit within the byte 'index'.

  if (nbits <= 8) { // Add the low-end bits to the telemetry packet.
    // Clear all bits except for the significant ones.
    byte = data[0] & (0xFF>>(8-nbits));
    // Determine the position where the bit has to be inserted.
    index   = packet->current_bit / 8;
    modulus = packet->current_bit % 8;
    // Insert the bits to the internal buffer.
    if (modulus+nbits<=8) {
      packet->data[index] = packet->data[index] |
	(unsigned char)(byte<<(8-nbits-modulus));
    } else {
      packet->data[index] = packet->data[index] |
	(unsigned char)(byte>>(nbits+modulus-8));
      packet->data[index+1] = packet->data[index+1] |
	(unsigned char)(byte<<(16-nbits-modulus));
    }
    // Update the counter for the current bit in the TelemetryPacket.
    packet->current_bit += nbits;

  } else { // The requested number of bits is bigger than 8.
           // Add the first n bytes to the Telemetry Packet and the high-
           // end bits from the last byte.
    for (index=0; index<nbits/8; index++) {
      addData2TelemetryPacket(packet, &(data[index]), 8);
    }
    modulus = nbits % 8;
    if (0!=modulus) {
      byte = data[index]>>(8-modulus);
      addData2TelemetryPacket(packet, &byte, modulus);
    }
  }

  return(EXIT_SUCCESS);
}



int writeTelemetryPacket2File(TelemetryPacket* packet, FILE* file)
{
  // Write the bytes from the internal storage to the file.
  int nbytes = fwrite (packet->data, 1, packet->nbits/8, file);
  if (nbytes < packet->nbits/8) {
    HD_ERROR_THROW("Error: Writing TelemetryPacket to binary output file failed!\n",
		   EXIT_FAILURE);
    return(EXIT_FAILURE);
  } else {
    return(EXIT_SUCCESS);
  }
}
