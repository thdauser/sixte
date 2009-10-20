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


/** Initialize the TelemetryPacket data structure.
 * The routine allocates the necessary memory and sets the values of 
 * the required properties. 
 * The return value of the function is the error status. */
int initTelemetryPacket(TelemetryPacket* packet, int nbits);

/** Release the memory allocated in the initTelemetryPacket() routine. */
void cleanupTelemetryPacket(TelemetryPacket* packet);

/** Clean the data in the TelemetryPacket internal storage 
 * (set all values to zero) and reset the internat bit counter. */
void newTelemetryPacket(TelemetryPacket* packet);

/** Returns the number of still unused and therefore available bits in the
 * TelemetryPacket. */
int availableBitsInTelemetryPacket(TelemetryPacket* packet);

/** Add binary data to the TelemetryPacket. 
 * The amount of 'nbits' is written from the memory section 'data' is pointing to,
 * and is inserted at the position designated by 'current_bit' in the 
 * TelemetryPackets internal data storage. 
 * The return value of the function is the error status. */
int addData2TelemetryPacket(TelemetryPacket* packet, unsigned char* data, int nbits);


#endif /* TELEMETRYPACKET_H */
