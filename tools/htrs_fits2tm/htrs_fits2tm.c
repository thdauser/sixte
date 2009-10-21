#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include "sixt.h"
#include "htrsevent.h"
#include "htrseventfile.h"
#include "telemetrypacket.h"

#define TOOLSUB htrs_fits2tm_main
#include "headas_main.c"


// Program parameters.
struct Parameters {
  char eventlist_filename[FILENAME_LENGTH];
  char output_filename[FILENAME_LENGTH];

  // Number of bits in each telemetry packet.
  int n_packet_bits;
  // Number of bits for the telemetry header.
  int n_header_bits;
  // Number of bits per spectral bin.
  int n_bin_bits;
  // Number of response channels binned to one spectral bin.
  int bin_size;
  // Integration time for one spectrum.
  double integration_time;
};


// Reads the program parameters using PIL
static int getpar(struct Parameters* parameters);


//////////////////////////////////
//    MAIN
int htrs_fits2tm_main()
{
  // Number of channels in the response file.
  const long n_rsp_channels = 1024;
  // Program parameters.
  struct Parameters parameters;

  HTRSEventFile eventfile;
  FILE* output_file=NULL;
  TelemetryPacket tmpacket;

  // Number of bins in the spectrum.
  int n_spectrum_bins=0;
  // Buffer for the spectral binning.
  int* spectrum=NULL;
  // Number of bits per binned spectrum.
  int n_spectrum_bits=0;
  // Byte buffer for the spectrum. In order to write it to the TelemetryPacket
  // the spectrum is converted to a binary byte buffer.
  unsigned char* byte_buffer=NULL;
  // Due to the limited number of bits per spectral bin, the maximum number
  // of counts per bin is limited. The following variable represents this maximum.
  int max_counts=0;

  int count; // Counter.

  char msg[MAXMSG]; // Error message buffer.
  int status=EXIT_SUCCESS;


  // HEATOOLs: register program
  set_toolname("htrs_fits2tm");
  set_toolversion("0.01");


  do { // Beginning of ERROR handling loop

    // --- Initialization ---
 
    // Read parameters using PIL library:
    if ((status=getpar(&parameters))) break;

    // Open the event list FITS file.
    status=openHTRSEventFile(&eventfile, parameters.eventlist_filename, READONLY);
    if(EXIT_SUCCESS!=status) break;

    // Determine the number of bits required for one spectrum.
    if (n_rsp_channels % parameters.bin_size != 0) {
      n_spectrum_bins = n_rsp_channels/parameters.bin_size +1;
    } else {
      n_spectrum_bins = n_rsp_channels/parameters.bin_size;
    }
    n_spectrum_bits = n_spectrum_bins * parameters.n_bin_bits;
    // Currently the number of bits must be a multiple of 8, i.e., we are only
    // dealing with complete bytes.
    assert(n_spectrum_bits % 8 == 0);

    // Determine the maximum number of counts per spectral bin.
    unsigned char maximum=0;
    for (count=0; count<parameters.n_bin_bits; count++) {
      maximum = (maximum<<1) + 1;
    }
    max_counts = (int)maximum;
    printf("maximum number of counts per spectral bin: %d\n", max_counts);
    
    // Open the binary file for output:
    output_file = fopen(parameters.output_filename, "w+");
    if (NULL==output_file) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: output file '%s' could not be opened!\n", 
	      parameters.output_filename);
      HD_ERROR_THROW(msg, status);
      break;
    }

    // Allocate memory for the spectrum buffer.
    spectrum=(int*)malloc(n_spectrum_bins*sizeof(int));
    if (NULL==spectrum) {
      status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: Could not allocate memory for spectrum buffer!\n", status);
      break;
    }
    for(count=0; count<n_spectrum_bins; count++) {
      spectrum[count]=0;
    }

    // Allocate memory for the byte buffer.
    byte_buffer=(unsigned char*)malloc(n_spectrum_bits/8*sizeof(unsigned char*));
    if (NULL==byte_buffer) {
      status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: Could not allocate memory for byte buffer!\n", status);
      break;
    }

    // Initialize the TelemetryPacket data structure.
    status=initTelemetryPacket(&tmpacket, parameters.n_packet_bits);
    if(EXIT_SUCCESS!=status) break;

    // Start a new (the first) Telemetry Packet.
    newTelemetryPacket(&tmpacket);
    // Add the Packet header.
    for (count=0; count*8<parameters.n_header_bits; count++) {
      byte_buffer[count] = 0;
    }
    status = addData2TelemetryPacket(&tmpacket, byte_buffer, count*8);
    if (EXIT_SUCCESS!=status) break;
    

    // --- END of Initialization ---


    // Loop repeated as long as there are entries in the event file.
    HTRSEvent event = { .time=0. };
    double spectrum_time=0.;
    while ((EXIT_SUCCESS==status)&&(0==EventFileEOF(&eventfile.generic))) {
      
      // Read the next event from the event file.
      status=HTRSEventFile_getNextRow(&eventfile, &event);

      // If the time exceed since starting the spectral binning is larger 
      // than the binning time, close the spectrum and add it to the 
      // current TelemetryPacket.
      if (event.time-spectrum_time > parameters.integration_time) {
	
	// Check whether the TelemetryPacket can store another spectrum,
	// or whether it is already full. In the latter case, write
	// the data to the binary output file and start a new packet.
	if (availableBitsInTelemetryPacket(&tmpacket)<n_spectrum_bits) {
	  // Store the telemetry packet in the output binary file.
	  writeTelemetryPacket2File(&tmpacket, output_file);

	  // Start a new TelemetryPacket.
	  newTelemetryPacket(&tmpacket);
	  // Add the Packet header.
	  for (count=0; count*8<parameters.n_header_bits; count++) {
	    byte_buffer[count] = 0;
	  }
	  status = addData2TelemetryPacket(&tmpacket, byte_buffer, count*8);
	  if (EXIT_SUCCESS!=status) break;
	}
	
	// Convert the spectrum to binary format and add it to the 
	// TelemetryPacket.
	int byte_index, bit_in_byte;
	for(byte_index=0; byte_index<(n_spectrum_bits/8); byte_index++) {
	  byte_buffer[byte_index]=0;
	}
	for(count=0; count<n_spectrum_bins; count++) {
	  // Check whether the maximum number of counts in this spectral bin is exceeded.
	  // (Overflows due to the limited number of bits per bin have to be avoided.)
	  if(spectrum[count]>max_counts) {
	    printf("Warning: overflow (maximum number of counts per bin exceeded)!\n");
	    spectrum[count]=max_counts;
	  }

	  byte_index = (count*parameters.n_bin_bits)/8;
	  bit_in_byte = (count*parameters.n_bin_bits)%8;
	  // Check whether the current spectral bin fits within the current
	  // byte, ...
	  if (bit_in_byte <= 8-parameters.n_bin_bits) {
	    byte_buffer[byte_index] += 
	      (unsigned char)(spectrum[count]<<(8-parameters.n_bin_bits-bit_in_byte));
	  } else {
	    // ... or whether it overlaps with 2 bytes.
	    byte_buffer[byte_index] += 
	      (unsigned char)(spectrum[count]>>(parameters.n_bin_bits+bit_in_byte-8));
	    byte_buffer[byte_index+1] += 
	      (unsigned char)(spectrum[count]<<(16-parameters.n_bin_bits-bit_in_byte));
	  }

	  // Clear the spectrum buffer.
	  spectrum[count]=0;
	}
	
	status = addData2TelemetryPacket(&tmpacket, byte_buffer, n_spectrum_bits);
	if(EXIT_SUCCESS!=status) break;

	// Update the time for the spectral binning.
	spectrum_time += parameters.integration_time;

      } else { 
	// The newly read event has to be added to the current spectrum.
	// (The binning time was not exceeded.)
	count = event.pha/parameters.bin_size;
	spectrum[count]++;

      } // END Newly read event belongs to current spectrum.
    } // END of loop over all entries in the event file.

  } while (0); // END of ERROR handling loop


  // --- Clean up ---

  // Release memory for the spectrum buffer.
  if (NULL!=spectrum) free(spectrum);

  // Release the memory for the internal storage of the TelemetryPacket.
  cleanupTelemetryPacket(&tmpacket);

  // Close files.
  if (NULL!=output_file) fclose(output_file);
  closeHTRSEventFile(&eventfile);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully\n\n");
  return(status);
}



static int getpar(struct Parameters* parameters)
{
  int status=EXIT_SUCCESS;

  // Get the name of the input file (FITS event list).
  if ((status = PILGetFname("eventlist_filename", parameters->eventlist_filename))) {
    HD_ERROR_THROW("Error reading the filename of the input event file (FITS)!\n", status);
  }

  // Get the name of the output file (binary).
  else if ((status = PILGetFname("output_filename", parameters->output_filename))) {
    HD_ERROR_THROW("Error reading the filename of the output file (binary)!\n", status);
  }

  else if ((status = PILGetInt("n_packet_bits", &parameters->n_packet_bits))) {
    HD_ERROR_THROW("Error reading the number of bits per telemetry packet!\n", status);
  }

  else if ((status = PILGetInt("n_header_bits", &parameters->n_header_bits))) {
    HD_ERROR_THROW("Error reading the number of bits for the packet header!\n", status);
  }

  else if ((status = PILGetInt("n_bin_bits", &parameters->n_bin_bits))) {
    HD_ERROR_THROW("Error reading the number of bits per spectral bin!\n", status);
  }

  else if ((status = PILGetInt("bin_size", &parameters->bin_size))) {
    HD_ERROR_THROW("Error reading the spectral bin size (response channels per bin)!\n", status);
  }

  else if ((status = PILGetReal("integration_time", &parameters->integration_time))) {
    HD_ERROR_THROW("Error reading the integration time per spectrum!\n", status);
  }

  return(status);
}

