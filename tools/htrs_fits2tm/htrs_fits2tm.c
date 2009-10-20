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
  FILE *output_file=NULL;
  TelemetryPacket tmpacket;

  //  unsigned int max = 0;

  char msg[MAXMSG];    // error message buffer
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

    // Open the binary file for output:
    output_file = fopen(parameters.output_filename, "w+");
    if (NULL==output_file) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: output file '%s' could not be opened!\n", 
	      parameters.output_filename);
      HD_ERROR_THROW(msg, status);
      break;
    }

    // Initialize the TelemetryPacket data structure.
    status=initTelemetryPacket(&tmpacket, parameters.n_packet_bits);
    if(EXIT_SUCCESS!=status) break;


      /*  } else if (mode == MODE_SPECTRUM) { 
    // The individual events are binned to spectra before they are converted
    // to a byte stream.

    struct eROSITAEvent event;
    event.pha = 0;
    event.xi = 0;
    event.yi = 0;
    event.patid = 0;
    event.patnum = 0;
    event.frame = 0;
    event.pileup = 0;

    double time = binning_time; // Time of the current spectrum (end of spectrum)
    long channel;      // Channel counter
    unsigned char spectrum[Nchannels];  // Buffer for binning the spectrum.
    unsigned char output_buffer[N_HTRS_BYTES];
    
    // Clear the output buffer:
    binary_output_clear_bytes(output_buffer, N_HTRS_BYTES);
    // Clear the spectrum:
    binary_output_clear_bytes(spectrum, Nchannels);
    
    for (eventlist_file.row=0; eventlist_file.row<eventlist_file.nrows;
	 eventlist_file.row++) {
      
      // Read the event from the FITS file.
      if (get_eventlist_row(eventlist_file, &event, &status)) break;
      
      // Check whether binning time was exceeded:
      if (event.time > time) {
	// Store binned spectrum, clear spectrum buffer and start
	// new binning cycle:
	for (channel=0; channel<Nchannels; channel+=2) {	  
	  
	  if (spectrum[channel]   > max) max = spectrum[channel];
	  if (spectrum[channel+1] > max) max = spectrum[channel+1];
	  
	  output_buffer[9 + (channel/2)%119] = (unsigned char)
	    ((spectrum[channel] << 4) + (spectrum[channel+1] & 0x0F));
	  
	  // Clear binned spectrum:
	  spectrum[channel]   = 0;  
	  spectrum[channel+1] = 0;
	  
	  if ((channel+2)%28 == 0) { // Byte frame (128 byte) is complete!
	    
	    // Syncword 1 and 2:
	    // output_buffer[0] = (char)'K';
	    // output_buffer[1] = (char)'R';
	    output_buffer[0] = 0x4B;  // 'K'
	    output_buffer[1] = 0x82;  // 'R'
	    
	      // Spectrum Time:
	    long ltime = (long)(time/binning_time);
	    output_buffer[2] = (unsigned char)(ltime>>24);
	    output_buffer[3] = (unsigned char)(ltime>>16);
	    output_buffer[4] = (unsigned char)(ltime>>8);
	    output_buffer[5] = (unsigned char)ltime;
	    //headas_chat(5, "%ld: %u %u %u %u\n", ltime, output_buffer[2], 
	    //	output_buffer[3], output_buffer[4], output_buffer[5]);
	    
	    // Spectrum Sequence counter:
	    output_buffer[6] = (channel/238);

	    // Data type ID:
	    output_buffer[7] = 0x83;  // 'S'
	    // 0x..   -> hexadecimal
	    // 0...   -> octal
	    
	    // Number of used bytes:
	    if (channel/119 == 4) {
	      output_buffer[8] = 0x24; // 36
	    } else {
	      output_buffer[8] = 0x77; // 119
	    }
	    
	    // Write bytes to file
	    int nbytes = fwrite (output_buffer, 1, N_HTRS_BYTES, output_file);
	    if (nbytes < N_HTRS_BYTES) {
	      status=EXIT_FAILURE;
	      sprintf(msg, "Error: writing data to output file '%s' failed!\n", 
		      output_filename);
	      HD_ERROR_THROW(msg,status);
	    }
	    
	    // Clear the output buffer:
	    binary_output_clear_bytes(output_buffer, N_HTRS_BYTES);
	    
	  } // END of starting new byte frame
	  
	} // END of loop over binned spectrum
	time += binning_time;  // next binning cycle

      }

      if ((event.pha<=0)||(event.pha>Nchannels)) printf("Error!!\n");
      
      // Add event to the spectrum.
      spectrum[event.pha-1]++;
      
    } // END of loop over all events in the FITS file
    
    headas_chat(5, "maximum spectral bin: %u\n", max);
  
    } // END (mode == MODE_SPECTRUM) */

  } while (0); // END of ERROR handling loop


  // --- Clean up ---

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

