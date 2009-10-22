#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include "sixt.h"
#include "htrsevent.h"
#include "htrseventfile.h"
#include "htrstelemetrypacket.h"

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
  HTRSTelemetryPacket htmpacket;

  // Number of bins in the spectrum.
  int n_spectrum_bins=0;
  // Buffer for the spectral binning.
  int* spectrum=NULL;

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

    // Determine the number of spectral bins.
    if (n_rsp_channels % parameters.bin_size != 0) {
      n_spectrum_bins = n_rsp_channels/parameters.bin_size +1;
    } else {
      n_spectrum_bins = n_rsp_channels/parameters.bin_size;
    }
    
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

    // Initialize the TelemetryPacket data structure.
    struct HTRSTelemetryPacketParameters htpp = {
      .n_packet_bits = parameters.n_packet_bits,
      .n_header_bits = parameters.n_header_bits,
      .n_spectrum_bits = n_spectrum_bins * parameters.n_bin_bits,
      .n_bin_bits      = parameters.n_bin_bits
    };
    status=initHTRSTelemetryPacket(&htmpacket, &htpp);
    if(EXIT_SUCCESS!=status) break;

    // Start a new (the first) Telemetry Packet.
    status = newHTRSTelemetryPacket(&htmpacket);
    if (EXIT_SUCCESS!=status) break;
    

    // --- END of Initialization ---


    // Loop repeated as long as there are entries in the event file.
    HTRSEvent event = { .time=0. };
    double spectrum_time=0.;
    while ((EXIT_SUCCESS==status)&&(0==EventFileEOF(&eventfile.generic))) {
      
      // Read the next event from the event file.
      status=HTRSEventFile_getNextRow(&eventfile, &event);

      // If the time exceeded since starting the spectral binning is larger 
      // than the binning time, close the spectrum and add it to the 
      // current TelemetryPacket.
      if (event.time-spectrum_time > parameters.integration_time) {
	
	// Check whether the TelemetryPacket can store another spectrum,
	// or whether it is already full. In the latter case, write
	// the data to the binary output file and start a new packet.
	if (availableBitsInHTRSTelemetryPacket(&htmpacket)<htmpacket.n_spectrum_bits) {
	  // Store the telemetry packet in the output binary file.
	  writeHTRSTelemetryPacket2File(&htmpacket, output_file);

	  // Start a new TelemetryPacket.
	  status = newHTRSTelemetryPacket(&htmpacket);
	  if (EXIT_SUCCESS!=status) break;
	}
	
	// Convert the spectrum to binary format and add it to the 
	// TelemetryPacket.
	status = addSpectrum2HTRSTelemetryPacket(&htmpacket, spectrum, n_spectrum_bins);
	if(EXIT_SUCCESS!=status) break;

	// Clear the spectrum buffer.
	for(count=0; count<n_spectrum_bins; count++) {
	  spectrum[count]=0;
	}
	
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
  cleanupHTRSTelemetryPacket(&htmpacket);

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

