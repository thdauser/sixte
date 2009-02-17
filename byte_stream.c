#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>


#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#define TOOLSUB byte_stream_main
#include "headas_main.c"

#include "global_constants.h"
#include "event_list.h"



// This routine reads the program parameters using PIL.
int byte_stream_getpar(char input_filename[], char output_filename[],
		       double* binning_time);



//////////////////////////////////
int byte_stream_main()
{
  struct Event_List_File event_list_file;   // FITS file
  char output_filename[FILENAME_LENGTH];
  const long Nchannels = 1024;
  int spectrum[Nchannels];  // Buffer for binning the spectrum.
  double time = 0.0;   // Time of the current spectrum (end of spectrum)
  double binning_time; // Delta t (time step, length of each spectrum)

  char msg[MAXMSG];    // error message buffer
  int status=EXIT_SUCCESS;


  // HEATOOLs: register program
  set_toolname("byte_stream");
  set_toolversion("0.01");


  // Get the parameters:
  status = byte_stream_getpar(event_list_file.filename, output_filename, 
			      &binning_time);
  if (status != EXIT_SUCCESS) return (status);

  do { // Beginning of ERROR handling loop
    int hdutype;
    if (fits_open_table(&event_list_file.fptr, event_list_file.filename, 
			READONLY, &status)) break;

    // Get the HDU type.
    if (fits_get_hdu_type(event_list_file.fptr, &hdutype, &status)) break;
    
    // If the open HDU is an image extension, throw an error.
    if (hdutype==IMAGE_HDU) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: input file '%s' contains no binary table!\n",
	      event_list_file.filename);
      HD_ERROR_THROW(msg,status);
      break;
    }

    // Determine the number of rows in the event list.
    if (fits_get_num_rows(event_list_file.fptr, &event_list_file.nrows, &status)) 
      break;



    // Loop over all events in the FITS file:
    headas_chat(5, "processing events ...\n");
    struct Event event;
    long channel;      // channel counter

    time = binning_time;
    for (channel=0; channel<Nchannels; channel++) {
      spectrum[channel] = 0;
    }

    for (event_list_file.row=0; event_list_file.row<event_list_file.nrows;
	 event_list_file.row++) {
      
      if (get_eventtbl_row(event_list_file, &event, &status)) break;

      
      // Check whether binning time was exceeded:
      if (event.time > time) {
	// Store binned spectrum clear spectrum buffer and start 
	// new binning cycle:
	for (channel=0; channel<Nchannels; channel++) {
	  if (channel%119 == 0) { // new byte frame (128 byte)
	    fprintf(stdout, "\n");

	    // Syncword 1 and 2:
	    fprintf(stdout, "%c%c", 'K', 'R');
	    
	    // Spectrum Time:
	    fprintf(stdout, "%c%c%c%", 0, 0, 0, 0);

	    // Spectrum Sequence counter:
	    fprintf(stdout, "%c", channel/119);

	    // Data type ID:
	    fprintf(stdout, "%c", 'S');

	
	  }
	  
	  spectrum[channel] = 0;  // clear buffer
	}
	time += binning_time;  // next binning cycle
      }

      if ((event.pha<=0)||(event.pha>Nchannels)) printf("Error!!\n");

      // Add event to the spectrum.
      spectrum[event.pha-1]++;

    } // END of loop over all events in the FITS file

  } while (0); // END of ERROR handling loop


  return(status);
}





//////////////////////////////////////////////////////
int byte_stream_getpar(
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
