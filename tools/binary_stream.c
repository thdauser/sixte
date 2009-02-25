#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#define TOOLSUB binary_stream_main
#include "headas_main.c"

#include "sixt.h"
#include "detectors.h"



// This routine reads the program parameters using PIL.
int binary_stream_getpar(char input_filename[], char output_filename[],
		       double* binning_time);



//////////////////////////////////
int binary_stream_main()
{
  struct Eventlist_File eventlist_file;   // FITS file
  char output_filename[FILENAME_LENGTH];
  FILE* output_file = NULL;
  double binning_time; // Delta t (time step, length of each spectrum)

  const long Nchannels = 1024;
  unsigned int max = 0;

  char msg[MAXMSG];    // error message buffer
  int status=EXIT_SUCCESS;


  // HEATOOLs: register program
  set_toolname("binary_stream");
  set_toolversion("0.01");


  // Get the parameters:
  status = binary_stream_getpar(eventlist_file.filename, output_filename, 
			      &binning_time);
  if (status != EXIT_SUCCESS) return (status);

  do { // Beginning of ERROR handling loop
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
      HD_ERROR_THROW(msg,status);
      break;
    }

    // Determine the number of rows in the event list.
    if (fits_get_num_rows(eventlist_file.fptr, &eventlist_file.nrows, &status)) 
      break;


    // Open binary file for output:
    output_file = fopen(output_filename, "w+");
    if (output_file == NULL) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: output file '%s' could not be opened!\n", 
	      output_filename);
      HD_ERROR_THROW(msg,status);
      break;
    }




    // Loop over all events in the FITS file:
    headas_chat(5, "processing events ...\n");

    double time = 0.0;   // Time of the current spectrum (end of spectrum)

    struct Event event;
    event.pha = 0;
    event.xi = 0;
    event.yi = 0;
    event.grade = 0;
    event.patid = 0;
    event.patnum = 0;
    event.frame = 0;
    event.pileup = 0;

    long channel;    // Channel counter
    int count;       // Counter for access to output buffer
    unsigned int spectrum[Nchannels];  // Buffer for binning the spectrum.
    unsigned char output_buffer[128];

    // Clear the output buffer:
    for (count = 0; count < 128; count++) {
      output_buffer[count] = 0;
    }

    // Clear the spectrum:
    time = binning_time;
    for (channel=0; channel<Nchannels; channel++) {
      spectrum[channel] = 0;
    }

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

	  if ((channel+2)%238 == 0) { // Byte frame (128 byte) is complete!

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
	    int nbytes = fwrite (output_buffer, 1, 128, output_file);
	    if (nbytes < 128) {
	      status=EXIT_FAILURE;
	      sprintf(msg, "Error: writing data to output file '%s' failed!\n", 
		      output_filename);
	      HD_ERROR_THROW(msg,status);
	    }

	    // Clear the output buffer:
	    for (count = 0; count < 128; count++) {
	      output_buffer[count] = 0;
	    }
	    
	  } // END of starting new byte frame

	} // END of loop over binned spectrum
	time += binning_time;  // next binning cycle

      }

      if ((event.pha<=0)||(event.pha>Nchannels)) printf("Error!!\n");

      // Add event to the spectrum.
      spectrum[event.pha-1]++;

    } // END of loop over all events in the FITS file

  } while (0); // END of ERROR handling loop


  headas_chat(5, "maximum spectral bin: %u\n", max);


  // Close files
  if (output_file) fclose(output_file);
  if (eventlist_file.fptr) fits_close_file(eventlist_file.fptr, &status);

  return(status);
}





//////////////////////////////////////////////////////
int binary_stream_getpar(
			 char input_filename[],  // FITS event list
			 char output_filename[], // binary output file
			 double* binning_time    // time span for spectral binning
			 )
{
  int status=0;        // error status
  char msg[MAXMSG];    // error message buffer

  if ((status = PILGetFname("input_filename", input_filename))) {
    sprintf(msg, "Error reading the filename of the input event list (FITS)!\n");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetFname("output_filename", output_filename))) {
    sprintf(msg, "Error reading the filename of the output file (binary)!\n");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetReal("binning_time", binning_time))) {
    sprintf(msg, "Error reading the spectral binning time!\n");
    HD_ERROR_THROW(msg,status);
  }


  return(status);
}
