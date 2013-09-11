#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include "sixt.h"
#include "htrsevent.h"
#include "htrseventfile.h"
#include "htrstelstream.h"

#define TOOLSUB htrs_fits2tm_main
#include "headas_main.c"


// Program parameters.
struct Parameters {
  char eventlist_filename[MAXFILENAME];
  char output_filename[MAXFILENAME];

  // Number of bits in each telemetry packet.
  int n_packet_bits;
  // Number of bits for the telemetry header.
  int n_header_bits;
  // Number of bits per spectral bin.
  int n_bin_bits;
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
  HTRSTelStream* htstream=NULL;

  int status=EXIT_SUCCESS;


  // HEATOOLs: register program
  set_toolname("htrs_fits2tm");
  set_toolversion("0.01");


  do { // Beginning of ERROR handling loop

    // --- Initialization ---
 
    // Read parameters using PIL library:
    if ((status=getpar(&parameters))) break;

    // Open the event list FITS file.
    status=openHTRSEventFile(&eventfile, parameters.eventlist_filename, 
			     READONLY);
    if(EXIT_SUCCESS!=status) break;

    // Set up the look-up table that assigns a spectrum bin to each 
    // detector RSP channel.
    int chans2bins[n_rsp_channels];
    int count, count2;
    for (count=0; count<128; count++) {
      chans2bins[count] = count;
      for (count2=0; count2<7; count2++) {
	chans2bins[128+count*7+count2] = count+128;
      }
    }

    // Initialize the TelemetryPacket data structure.
    struct HTRSTelStreamParameters htsp = {
      .n_packet_bits = parameters.n_packet_bits,
      .n_header_bits = parameters.n_header_bits,
      .n_bin_bits    = parameters.n_bin_bits,
      .n_channels = n_rsp_channels,
      .n_bins     = 256,
      .integration_time = parameters.integration_time,
      .chans2bins = chans2bins,
      .output_filename = parameters.output_filename // Copy string address !!
    };
    htstream=getHTRSTelStream(&htsp, &status);
    if(EXIT_SUCCESS!=status) break;

    // --- END of Initialization ---


    // Loop repeated as long as there are entries in the event file.
    HTRSEvent event = { .time=0. };
    while ((EXIT_SUCCESS==status)&&(0==EventListEOF(&eventfile.generic))) {
      
      // Read the next event from the event file.
      status=HTRSEventFile_getNextRow(&eventfile, &event);
      if (EXIT_SUCCESS!=status) break;

      // Add the event to the HTRSTelStream.
      status=addEvent2HTRSTelStream(htstream, &event);
      if (EXIT_SUCCESS!=status) break;

    } // END of loop over all entries in the event file.
    if (EXIT_SUCCESS!=status) break;
    
    // Finalize the HTRSTelStream. Otherwise the current TelemetryPacket will
    // not be written to the output file and will therefore be lost.
    status = finalizeHTRSTelStream(htstream);

    // Print the statistical information about the generated Telemetry Stream.
    HTRSTelStreamPrintStatistics(htstream);

  } while (0); // END of ERROR handling loop


  // --- Clean up ---

  // Release the memory for the internal storage of the TelemetryPacket.
  freeHTRSTelStream(htstream);

  // Close files.
  closeHTRSEventFile(&eventfile);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully\n\n");
  return(status);
}



static int getpar(struct Parameters* parameters)
{
  int status=EXIT_SUCCESS;

  // Get the name of the input file (FITS event list).
  if ((status = PILGetFname("eventlist_filename", 
			    parameters->eventlist_filename))) {
    HD_ERROR_THROW("Error reading the filename of the input event file "
		   "(FITS)!\n", status);
  }

  // Get the name of the output file (binary).
  else if ((status = PILGetFname("output_filename", 
				 parameters->output_filename))) {
    HD_ERROR_THROW("Error reading the filename of the output file "
		   "(binary)!\n", status);
  }

  else if ((status = PILGetInt("n_packet_bits", &parameters->n_packet_bits))) {
    HD_ERROR_THROW("Error reading the number of bits per telemetry packet!\n", 
		   status);
  }

  else if ((status = PILGetInt("n_header_bits", &parameters->n_header_bits))) {
    HD_ERROR_THROW("Error reading the number of bits for the packet header!\n",
		   status);
  }

  else if ((status = PILGetInt("n_bin_bits", &parameters->n_bin_bits))) {
    HD_ERROR_THROW("Error reading the number of bits per spectral bin!\n", 
		   status);
  }

  else if ((status = PILGetReal("integration_time", 
				&parameters->integration_time))) {
    HD_ERROR_THROW("Error reading the integration time per spectrum!\n", 
		   status);
  }

  return(status);
}

