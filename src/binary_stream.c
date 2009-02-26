#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <malloc.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#define TOOLSUB binary_stream_main
#include "headas_main.c"

#include "sixt.h"
#include "detectors.h"


// Number of bytes per TLM byte frame
#define N_HTRS_BYTES    (128)
#define N_EROSITA_BYTES (128)

#include "byte_output.c"


///////////////////////////////////////////////////
// Convert a squence of chars into captial letters. 
// The squence has to be terminated by a '\0' mark.
void strtoupper(char string[]) {
  int count=0;
  while (string[count] != '\0') {
    string[count] = toupper(string[count]);
  };
}




//////////////////////////////////
//    MAIN
int binary_stream_main()
{
  // Program output mode (events or spectrum)
  enum Mode {
    MODE_INVALID =0,
    MODE_EVENTS  =1,
    MODE_SPECTRUM=2
  };
  enum Mode mode;

  struct Eventlist_File eventlist_file;   // FITS file
  char output_filename[FILENAME_LENGTH];
  FILE *output_file = NULL;
  double binning_time; // Delta t (time step, length of each spectrum)

  const long Nchannels = 1024;
  unsigned int max = 0;

  char msg[MAXMSG];    // error message buffer
  int status=EXIT_SUCCESS;


  // HEATOOLs: register program
  set_toolname("binary_stream");
  set_toolversion("0.01");


  do { // Beginning of ERROR handling loop

    // --- Initialization ---

    if ((status = PILGetFname("eventlist_filename", eventlist_file.filename))) {
      sprintf(msg, "Error reading the name of the input event list file (FITS)!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }

    int hdutype;
    eventlist_file.fptr=NULL;
    if (fits_open_table(&eventlist_file.fptr, eventlist_file.filename, 
			READONLY, &status)) break;

    // Get the HDU type.
    if (fits_get_hdu_type(eventlist_file.fptr, &hdutype, &status)) break;
    // If the open HDU is an image extension, throw an error.
    if (hdutype==IMAGE_HDU) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: input file '%s' contains no binary table!\n",
	      eventlist_file.filename);
      HD_ERROR_THROW(msg, status);
      break;
    }

    // Determine the output mode (events or spectrum) according to the 
    // telescope and detector type specified in the FITS header keywords.
    char telescop[MAXMSG], instrume[MAXMSG];
    char comment[MAXMSG]; // buffer
    if (fits_read_key(eventlist_file.fptr, TSTRING, "TELESCOP", telescop, 
		      comment, &status)) break;
    if (fits_read_key(eventlist_file.fptr, TSTRING, "INSTRUME", instrume, 
		      comment, &status)) break;

    // convert to captial letters:
    strtoupper(telescop); 
    strtoupper(instrume); 

    headas_chat(5, "TELESCOP: %s\n INSTRUME: %s\n", telescop, instrume);
    if (strcmp(telescop, "EROSITA") == 0) {
      headas_chat(5, "MODE: events\n");
      mode = MODE_EVENTS;
    } else if ((strcmp(telescop, "IXO") == 0) && 
	       (strcmp(instrume, "HTRS") == 0)) {
      headas_chat(5, "MODE: spectrum\n");
      mode = MODE_SPECTRUM;
    } else {
      mode = MODE_INVALID;
    }

    // Determine the number of rows in the event list.
    if (fits_get_num_rows(eventlist_file.fptr, &eventlist_file.nrows, &status)) 
      break;


    // Get the name of the output file (binary).
    if ((status = PILGetFname("output_filename", output_filename))) {
      sprintf(msg, "Error reading the filename of the output file (binary)!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }

    // If the output should be a spectrum determine the spectral binning time.
    if (mode == MODE_SPECTRUM) {
      if ((status = PILGetReal("binning_time", &binning_time))) {
	sprintf(msg, "Error reading the spectral binning time!\n");
	HD_ERROR_THROW(msg,status);
	break;
      }
    }

    // Open the binary file for output:
    output_file = fopen(output_filename, "w+");
    if (output_file == NULL) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: output file '%s' could not be opened!\n", 
	      output_filename);
      HD_ERROR_THROW(msg,status);
      break;
    }

    // --- END of Initialization ---


    // --- Beginning of EVENT PROCESSING ---

    // Loop over all events in the FITS file:
    headas_chat(5, "processing events ...\n");

    if (mode == MODE_EVENTS) {
      // EVENT mode, i.e., the events are transferred to a particular binary
      // data format without spectral binning or other modifications.
      
      struct Byte_Output *byte_output = get_Byte_Output(N_EROSITA_BYTES, output_file);
      struct Event *eventlist = (struct Event *)malloc(384*348*sizeof(struct Event));
      if ((byte_output == NULL)||(eventlist==NULL)) {
	status=EXIT_FAILURE;
	sprintf(msg, "Error: memory allocation failed!\n");
	HD_ERROR_THROW(msg,status);
	break;
      }
      int n_buffered_events=0;

      // Loop over all entries in the event list:
      for (eventlist_file.row=0; eventlist_file.row<eventlist_file.nrows;
	   eventlist_file.row++) {
      
	// Read the event from the FITS file.
	if (get_eventlist_row(eventlist_file, &(eventlist[n_buffered_events]), 
			      &status)) break;

	if (eventlist[n_buffered_events].frame > eventlist[0].frame) {
	  // Write the events to the binary output.
	  int count;
	  for (count=0; count<n_buffered_events; count++) {
	    if (byte_output_erosita_insert_event(byte_output, &(eventlist[count]))) {
	      status=EXIT_FAILURE;
	      sprintf(msg, "Error: generation of binary format failed!\n");
	      HD_ERROR_THROW(msg,status);
	      break;
	    }
	  }

	  if (byte_output_erosita_finish_frame(byte_output, eventlist[0].time)) {
	    status=EXIT_FAILURE;
	    sprintf(msg, "Error: generation of binary format failed!\n");
	    HD_ERROR_THROW(msg,status);
	    break;
	  }
	  
	  // New buffering period has started.
	  eventlist[0] = eventlist[n_buffered_events];
	  n_buffered_events = 0;

	} // END of loop over all buffered events

	n_buffered_events++;
	
      } // END of loop over all entries in the event list.

      if (status == EXIT_SUCCESS) {
	if (byte_output_erosita_finish_frame(byte_output, eventlist[0].time)) {
	  status=EXIT_FAILURE;
	  sprintf(msg, "Error: generation of binary format failed!\n");
	  HD_ERROR_THROW(msg,status);
	  break;
	}
      }

      if (status != EXIT_SUCCESS) {
	free_Byte_Output(byte_output);
	if (eventlist) free (eventlist);
	break;
      }

    } else if (mode == MODE_SPECTRUM) {
      // The individual events are binned to spectra before they are converted
      // to a byte stream.

      struct Event event;
      event.pha = 0;
      event.xi = 0;
      event.yi = 0;
      event.grade = 0;
      event.patid = 0;
      event.patnum = 0;
      event.frame = 0;
      event.pileup = 0;

      double time = binning_time; // Time of the current spectrum (end of spectrum)
      long channel;      // Channel counter
      unsigned char spectrum[Nchannels];  // Buffer for binning the spectrum.
      unsigned char output_buffer[N_HTRS_BYTES];

      // Clear the output buffer:
      byte_output_clear_bytes(output_buffer, N_HTRS_BYTES);
      // Clear the spectrum:
      byte_output_clear_bytes(spectrum, Nchannels);

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
	      byte_output_clear_bytes(output_buffer, N_HTRS_BYTES);
	    
	    } // END of starting new byte frame

	  } // END of loop over binned spectrum
	  time += binning_time;  // next binning cycle

	}

	if ((event.pha<=0)||(event.pha>Nchannels)) printf("Error!!\n");

	// Add event to the spectrum.
	spectrum[event.pha-1]++;

      } // END of loop over all events in the FITS file

      headas_chat(5, "maximum spectral bin: %u\n", max);
  
    } // END (mode == MODE_SPECTRUM)

  } while (0); // END of ERROR handling loop



  // Close files
  if (output_file) fclose(output_file);
  if (eventlist_file.fptr) fits_close_file(eventlist_file.fptr, &status);

  return(status);
}



