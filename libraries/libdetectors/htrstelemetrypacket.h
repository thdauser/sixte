#ifndef HTRSTELEMETRYPACKET_H
#define HTRSTELEMETRYPACKET_H 1

#include "sixt.h"
#include "telemetrypacket.h"


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


typedef struct {

  /** TelemetryPacket containing the actual data. */
  TelemetryPacket* tp;

  /** Number of bits reserved for the packet header.
   * This number is included in value nbits. */
  int n_header_bits;
  /** Number of bits per binned spectrum. */
  int n_spectrum_bits;
  /** Number of bits per spectrum bin. */
  int n_bin_bits;

  /** Maximum number of counts per spectral bin.
   * Due to the limited number of bits per spectral bin, the maximum number
   * of counts per bin is limited. The following variable represents this maximum. */
  int max_counts;

  /** Buffer for conversions to binary data. */
  unsigned char* byte_buffer;

} HTRSTelemetryPacket;


struct HTRSTelemetryPacketParameters {
  int n_packet_bits;
  int n_header_bits;
  int n_spectrum_bits;
  int n_bin_bits;
};


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////


/** Initialize the HTRSTelemetryPacket data structure.
 * The routine allocates the necessary memory and sets the values of 
 * the required properties. 
 * The return value of the function is the error status. */
int initHTRSTelemetryPacket(HTRSTelemetryPacket* packet, 
			    struct HTRSTelemetryPacketParameters* parameters);

/** Release the memory allocated in the initHTRSTelemetryPacket() routine. */
void cleanupHTRSTelemetryPacket(HTRSTelemetryPacket* packet);

/** Clean the data in the HTRSTelemetryPacket internal storage 
 * (set all values to zero) and reset the internat bit counter. 
 * The return value is the error status. */
int newHTRSTelemetryPacket(HTRSTelemetryPacket* packet);

/** Returns the number of still unused and therefore available bits in the
 * HTRSTelemetryPacket. */
int availableBitsInHTRSTelemetryPacket(HTRSTelemetryPacket* packet);

/** Add binary data to the HTRSTelemetryPacket. 
 * The amount of 'nbits' is written from the memory section 'data' is pointing to,
 * and is inserted at the position designated by 'current_bit' in the 
 * HTRSTelemetryPackets internal data storage. 
 * The return value of the function is the error status. */
int addSpectrum2HTRSTelemetryPacket(HTRSTelemetryPacket* packet, int* spectrum, 
				    int n_spectrum_bins);

/** Write the data stored in the HTRSTelemetryPacket to an output stream 
 * (usually a file).
 * The return value of the function is the error status. */
int writeHTRSTelemetryPacket2File(HTRSTelemetryPacket* packet, FILE* file);


#endif /* HTRSTELEMETRYPACKET_H */

