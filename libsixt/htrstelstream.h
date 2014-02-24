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
*/

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
  /** Number of bits per spectrum bin. */
  int n_bin_bits;

  /** Number of RSP Channels. */
  int n_channels; 
  /** Number of bins for the spectrum. */
  int n_bins; 

  /** Maximum number of counts per spectral bin. Due to the limited
   * number of bits per spectral bin, the maximum number of counts per
   * bin is limited. The following variable represents this
   * maximum. */
  int max_counts;
  /** Number of bit overflows. */
  int n_overflows;

  /** Binary output file. The data added to the HTRSTelemetryStream is
   * written to this binary output file. */
  FILE* output_file;

  /** Buffer for the binned spectrum. The number of elements is given by
   * n_bins. */
  int* spectrum;
  /** Start time of spectrum. This value is used to determine whether
   * events belong the current binned spectrum, or whether a new one
   * has to be started. */
  double spectrum_time;
  /** Integration time for the binned spectrum. If the time of an
   * event is bigger than the sum of spectrum_time and
   * integration_time, it does not belong to the current spectrum, but
   * a new one has to be started. */
  double integration_time;

  /** Lookup table, which RSP channels are binned to which spectrum bins. 
   * The array contains n_channels elements. */
  int* chans2bins; 
  
  /** Total number of spectra written to a telemetry packet. */
  int n_spectra;
  /** Total number of generated (maybe incomplete) packets. */
  int n_packets;
  /** Number of packets written to binary output file. */
  int n_written_packets;

} HTRSTelStream;


struct HTRSTelStreamParameters {
  int n_packet_bits;
  int n_header_bits;
  int n_bin_bits;

  int n_channels; 
  int n_bins; 

  double integration_time;

  /** Filename of the binary output file. */
  char* output_filename;

  /** Lookup table, which RSP channels are binned to which spectrum bins. */
  int* chans2bins; 
};


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////


/** Constructor of the HTRSTelStream object. Allocate memory
    and set up the initial configuration. */
HTRSTelStream* getHTRSTelStream
(struct HTRSTelStreamParameters* parameters, int* status);

/** Destructor of the HTRSTelStream object. Release the
    allocated memory. */
void freeHTRSTelStream(HTRSTelStream* stream);

/** Start a new TelemetryPacket and insert the packet header at the
    beginning. The return value is the error status. */
int newHTRSTelStreamPacket(HTRSTelStream* stream);

/** Add a new event to the HTRSTelStream. The function adds HTRSEvent
    objects to the HTRSTelStream. The return value is the error
    status. */
int addEvent2HTRSTelStream(HTRSTelStream* stream, HTRSEvent* event);

/** Finalize the HTRSTelStream. This routine add the current spectrum
    to the TelemetryPacket and writes the whole TelemetryPacket to the
    binary output file. It does NOT start a new TelemetryPacket. The
    return value is the error status. */
int finalizeHTRSTelStream(HTRSTelStream* stream);

/** Add the binned spectrum stored in the HTRSTelStream data structure
    to the Telemetry packet. After adding the spectrum to the
    TelemetryPacket the function clears the spectrum to start the next
    binning period. The return value is the error status. */
int HTRSTelStreamAddSpec2Packet(HTRSTelStream* stream);

/** Print statistical information about the generated Telemetry of the
    HTRS instrument. */
void HTRSTelStreamPrintStatistics(HTRSTelStream* stream);


#endif /* HTRSTELSTREAM_H */

