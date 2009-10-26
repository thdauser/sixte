#ifndef HTRSTELSTREAM_H
#define HTRSTELSTREAM_H 1

#include "sixt.h"
#include "htrsevent.h"
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

  /** Maximum number of counts per spectral bin. Due to the limited
   * number of bits per spectral bin, the maximum number of counts per
   * bin is limited. The following variable represents this
   * maximum. */
  int max_counts;

  /** Buffer for conversions to binary data. */
  unsigned char* byte_buffer;

} HTRSTelStream;


struct HTRSTelStreamParameters {
  int n_packet_bits;
  int n_header_bits;
  int n_spectrum_bits;
  int n_bin_bits;
};


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////


/** Constructor of the HTRSTelStream object. Allocate memory
 * and set up the initial configuration. */
HTRSTelStream* getHTRSTelStream
(struct HTRSTelStreamParameters* parameters, int* status);

/** Destructor of the HTRSTelStream object. Release the
 * allocated memory. */
void freeHTRSTelStream(HTRSTelStream* stream);

/** Add a new event to the HTRSTelStream. */
int addEvent2HTRSTelStream(HTRSTelStream* stream, HTRSEvent* event, 
			   int* status);


#endif /* HTRSTELSTREAM_H */

