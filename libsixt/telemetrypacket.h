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

#ifndef TELEMETRYPACKET_H
#define TELEMETRYPACKET_H 1

#include "sixt.h"


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


typedef struct {

  /** Number of bits per telemetry packet. */
  int nbits;

  /** Current bit in the telemetry packet.
   * If an IO-operation is performed this is the default starting bit
   * for the operation. */
  int current_bit;

  /** Data of the telemetry packet. */
  unsigned char* data;

} TelemetryPacket;


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////


/** Constructor for the TelemetryPacket data structure.  The routine
 * allocates the necessary memory and sets the values of the required
 * properties. */
TelemetryPacket* getTelemetryPacket(int nbits, int* status);

/** Destructor for the TelemetryPacket data structure. Release the
 * memory allocated in the initTelemetryPacket() routine. */
void freeTelemetryPacket(TelemetryPacket* packet);

/** Initialize a new TelemetryPacket. Clean the data in the
 * TelemetryPacket internal storage (set all values to zero) and reset
 * the internat bit counter. */
void newTelemetryPacket(TelemetryPacket* packet);

/** Check number of still available bits. Returns the number of still
 * unused and therefore available bits in the TelemetryPacket. */
int availableBitsInTelemetryPacket(TelemetryPacket* packet);

/** Add binary data to the TelemetryPacket. The number of 'nbits' bits
 * is of the input byte 'data' is inserted at the position designated
 * by 'current_bit' in the TelemetryPackets internal data storage.
 * <ul>
 * <li>If the number of bits is less or equal 8, the smallest
 * nbits of the byte data[0] are added to the packet. E.g., if
 * 'data[0]=0xF2' and 'nbits=3', the binary sequence '010' (=0x2) is
 * added to the TelemetryPacket.</li>
 * <li>If the number of bits is more than 8, the first bytes are
 * completely added to the TelemetryPacket, and the remaining bits are
 * taken from the beginning (high end) of the last byte. E.g., if
 * 'data[0]=0xA7', 'data[1]=0x5A', and 'nbits=12' the byte '0xA7' and
 * the bit sequence '0101' (=0x5) are added to the TelemetryPacket.</li>
 * </ul>
 * The return value of the
 * function is the error status. */
int addData2TelemetryPacket(TelemetryPacket* packet, unsigned char* data,
			    int nbits);

/** Data output to file stream. Write the data stored in the
 * TelemetryPacket to an output stream (usually a file). The return
 * value of the function is the error status. */
int writeTelemetryPacket2File(TelemetryPacket* packet, FILE* file);


#endif /* TELEMETRYPACKET_H */
