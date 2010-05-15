#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include "sixt.h"
#include "sixt_string.h"
#include "erositaevent.h"
#include "erositaeventfile.h"

#define TOOLSUB ero_fits2tm_main
#include "headas_main.c"


// Number of bytes per TLM byte frame
#define N_HTRS_BYTES    (128)
#define N_EROSITA_BYTES (128)


struct Parameters {
  char eventlist_filename[FILENAME_LENGTH];
};



// Data structure for the byte output writing routine.
struct Binary_Output {
  unsigned char *bytes;  // Byte output buffer
  FILE *fptr;            // Pointer to output file
  int n_bytes;           // Counting the bytes already written to the buffer
  int max_bytes;         // Maximum number of bytes in one TLM record
  unsigned char framecounter;  // Counter for the records (NOT detector frames)
  unsigned char framenumber;   // Counter for the records belonging to one and
                               // the same detector frame
};



// Clear the given byte buffer with 'length' bytes.
void binary_output_clear_bytes(unsigned char *bytes, const int length) {
  int count;
  for (count=0; count<length; count++) {
    bytes[count] = 0;
  }
}



// Constructor of Binary_Output
struct Binary_Output *get_Binary_Output(const int max_bytes, FILE *fptr) {
  // Get memory for the object itself:
  struct Binary_Output *binary_output = 
    (struct Binary_Output *) malloc(sizeof(struct Binary_Output));
  if (binary_output==NULL) return(NULL);

  // Get memory for the byte output buffer and clear it:
  binary_output->bytes = (unsigned char *) malloc(max_bytes*sizeof(unsigned char));
  if (binary_output->bytes==NULL) return(NULL);
  binary_output_clear_bytes(binary_output->bytes, max_bytes);

  // Set object configuration:
  binary_output->n_bytes=12; // -> byte address of the first element in the record
  binary_output->max_bytes=max_bytes;
  binary_output->framecounter=0;
  binary_output->framenumber=0;
  binary_output->fptr=fptr;

  return(binary_output);
}


// Destructor of Binary_Output
void free_Binary_Output(struct Binary_Output *binary_output) {
  // Free all allocated memory:
  if (binary_output) {
    if (binary_output->bytes) {
      free(binary_output->bytes);
    }
    free(binary_output);
  }
  binary_output=NULL;
}



// This routine finihes a record and writes all bytes to the output file.
int binary_output_erosita_finish_record(struct Binary_Output *binary_output){
  int count;

  // If the record already contains some elements but is not full,
  // fill it with 0x00 entries.
  if (binary_output->n_bytes>12) {
    for (count=binary_output->n_bytes; count<binary_output->max_bytes; count++) {
      binary_output->bytes[count] = 0x00;
    }

    // Write the header of the record:
    // syncwords
    binary_output->bytes[0] = 0xEB;
    binary_output->bytes[1] = 0x90;
    // frametime MSB
    binary_output->bytes[2] = 0x00;  // TODO
    binary_output->bytes[3] = 0x00;
    binary_output->bytes[4] = 0x00;
    binary_output->bytes[5] = 0x00;
    // framecounter
    binary_output->bytes[6] = binary_output->framecounter++;
    // header
    binary_output->bytes[7] = 0x1F;
    // DSP frame ID
    binary_output->bytes[8] = 0x20;
    // # of 16bit elements
    binary_output->bytes[9] = (unsigned char)((binary_output->n_bytes-12)/2);
    // # of total frame   ???????????  // TODO
    binary_output->bytes[10] = 0x00;
    // actual frame #
    binary_output->bytes[11] = binary_output->framenumber;


    // Write the byte buffer to the output file:
    int nbytes = fwrite (binary_output->bytes, 1, 
			 binary_output->max_bytes, binary_output->fptr);
    if (nbytes != binary_output->max_bytes) return(-1);
    

    // Reset variables:
    binary_output->n_bytes=12;

    // Clear output buffer:
    binary_output_clear_bytes(binary_output->bytes, binary_output->max_bytes);
  } // END (binary_output->n_bytes>12)

  return(0);
}



// Routine which is called to write an eROSITA event to the Binary_Output.
// Return value is '0' if everything is ok, otherwise the function returns '-1'.
int binary_output_erosita_insert_event(struct Binary_Output *binary_output,
				       eROSITAEvent *event)
{
  // Check for overflows:
  if (event->pha > 0x3FFF) return (-1);

  // Write the data of the event to the byte output buffer:
  binary_output->bytes[binary_output->n_bytes++] = 
    0x3F & (unsigned char)(event->pha>>8);
  binary_output->bytes[binary_output->n_bytes++] = 
    0xFF & (unsigned char)event->pha;
  binary_output->bytes[binary_output->n_bytes++] = 
    0xFF & (unsigned char)event->yi;
  binary_output->bytes[binary_output->n_bytes++] = 
    0xFF & (unsigned char)event->xi;

  // Check if record is full:
  if (binary_output->n_bytes>=binary_output->max_bytes) {
    if (binary_output_erosita_finish_record(binary_output)) return(-1);
  }

  return(0);
}



// This routine finihes a detector frame by adding the 3 last 4-byte elements
// containing the frame time, the number of discarded pixels and a spare element.
int binary_output_erosita_finish_frame(struct Binary_Output *binary_output, 
				     const double time) 
{
  // Write the TIME element (readout time of last detector frame):
  long ltime = (long)(time*40);
  if (ltime > 0x3FFFFFFF)  return(-1);
  binary_output->bytes[binary_output->n_bytes++] = 
    0x40 + (unsigned char)(ltime>>24);
  binary_output->bytes[binary_output->n_bytes++] = 
    (unsigned char)(ltime>>16);
  binary_output->bytes[binary_output->n_bytes++] = 
    (unsigned char)(ltime>>8);
  binary_output->bytes[binary_output->n_bytes++] = 
    (unsigned char)ltime;

  // Check if record is full:
  if (binary_output->n_bytes>=binary_output->max_bytes) {
    if (binary_output_erosita_finish_record(binary_output)) return(-1);
  }

  
  // Write number of DISCARDED PIXELS (assumed to be 0):
  binary_output->bytes[binary_output->n_bytes++] = 0xBE;
  binary_output->bytes[binary_output->n_bytes++] = 0x00;
  binary_output->bytes[binary_output->n_bytes++] = 0x00;
  binary_output->bytes[binary_output->n_bytes++] = 0x00;

  // Check if record is full:
  if (binary_output->n_bytes>=binary_output->max_bytes) {
    if (binary_output_erosita_finish_record(binary_output)) return(-1);
  }
  

  // Write SPARE element:
  binary_output->bytes[binary_output->n_bytes++] = 0xC0;
  binary_output->bytes[binary_output->n_bytes++] = 0x00;
  binary_output->bytes[binary_output->n_bytes++] = 0x00;
  binary_output->bytes[binary_output->n_bytes++] = 0x00;

  // Finish record in any case, even if it is not full yet:
  if (binary_output_erosita_finish_record(binary_output)) {
    return(-1);
  } else {
    binary_output->framenumber=0;
    return(0);
  }
}



//////////////////////////////////
//    MAIN
int ero_fits2tm_main()
{
  struct Parameters parameters;

  // Program output mode (events or spectrum)
  enum Mode {
    MODE_INVALID =0,
    MODE_EVENTS  =1,
    MODE_SPECTRUM=2
  };
  enum Mode mode;

  eROSITAEventFile eventlistfile; // FITS file containing the event list
  char output_filename[FILENAME_LENGTH];
  FILE *output_file = NULL;
  double binning_time; // Delta t (time step, length of each spectrum)

  char msg[MAXMSG];    // error message buffer
  int status=EXIT_SUCCESS;


  // HEATOOLs: register program
  set_toolname("ero_fits2tm");
  set_toolversion("0.01");


  do { // Beginning of ERROR handling loop

    // --- Initialization ---
 
    if ((status = PILGetFname("eventlist_filename", parameters.eventlist_filename))) {
      HD_ERROR_THROW("Error reading the name of the input event list file (FITS)!\n", status);
      break;
    }

    // Open the event list FITS file.
    status=openeROSITAEventFile(&eventlistfile, parameters.eventlist_filename, READONLY);
    if(EXIT_SUCCESS!=status) return(status);

    // Determine the output mode (events or spectrum) according to the 
    // telescope and detector type specified in the FITS header keywords.
    char telescop[MAXMSG], instrume[MAXMSG];
    char comment[MAXMSG]; // buffer
    if (fits_read_key(eventlistfile.generic.fptr, TSTRING, "TELESCOP", telescop, 
		      comment, &status)) break;
    if (fits_read_key(eventlistfile.generic.fptr, TSTRING, "INSTRUME", instrume, 
		      comment, &status)) break;

    // Convert to captial letters:
    strtoupper(telescop); 
    strtoupper(instrume); 

    headas_chat(5, "TELESCOP: %s\nINSTRUME: %s\n", telescop, instrume);
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

    
    /*if (mode == MODE_EVENTS) { */
    // EVENT mode, i.e., the events are transferred to a particular binary
    // data format without spectral binning or other modifications.
      
    struct Binary_Output *binary_output = 
      get_Binary_Output(N_EROSITA_BYTES, output_file);
    eROSITAEvent *eventlist = (eROSITAEvent*)malloc(10000*sizeof(eROSITAEvent));
    if ((NULL==binary_output)||(NULL==eventlist)) {
      status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation failed!\n", status);
      break;
    }
    // Clear the event buffer:
    int count;
    for (count=0; count<10000; count++) {
      eventlist[count].time=0.;
      eventlist[count].pha=0;
      eventlist[count].xi=0;
      eventlist[count].yi=0;
      eventlist[count].frame=0;
      eventlist[count].ra=SIXT_NAN;
      eventlist[count].dec=SIXT_NAN;
      eventlist[count].sky_xi=0;
      eventlist[count].sky_yi=0;	    
    }


    // Loop over all entries in the event list:
    int n_buffered_events=0;
    while(0==EventFileEOF(&eventlistfile.generic)) {

      // Read the event from the FITS file.
      status=eROSITAEventFile_getNextRow(&eventlistfile, &(eventlist[n_buffered_events]));
      if(EXIT_SUCCESS!=status) break;

      if (eventlist[n_buffered_events].frame > eventlist[0].frame) {
	// Write the events to the binary output.
	int count;
	for (count=0; count<n_buffered_events; count++) {
	  if (binary_output_erosita_insert_event(binary_output, &(eventlist[count]))) {
	    status=EXIT_FAILURE;
	    HD_ERROR_THROW("Error: generation of binary format failed!\n", status);
	    break;
	  }
	}
	// Finish the TLM record, even if it is not full yet.
	if (binary_output_erosita_finish_frame(binary_output, eventlist[0].time)) {
	  status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: generation of binary format failed!\n", status);
	  break;
	}
	  
	// New buffering period has started.
	eventlist[0] = eventlist[n_buffered_events];
	n_buffered_events = 0;

      } // END of loop over all buffered events

      n_buffered_events++;
	
    } // END of loop over all entries in the event list.

    if (EXIT_SUCCESS==status) {
      // Write the events to the binary output.
      int count;
      for (count=0; count<n_buffered_events; count++) {
	if (binary_output_erosita_insert_event(binary_output, &(eventlist[count]))) {
	  status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: generation of binary format failed!\n", status);
	  break;
	}
      }

      if (binary_output_erosita_finish_frame(binary_output, eventlist[0].time)) {
	status=EXIT_FAILURE;
	HD_ERROR_THROW("Error: generation of binary format failed!\n", status);
	break;
      }
    }
    if (EXIT_SUCCESS!=status) {
      free_Binary_Output(binary_output);
      if (eventlist) free(eventlist);
      break;
    }

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

  // Close files
  if (output_file) fclose(output_file);
  closeeROSITAEventFile(&eventlistfile);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully\n\n");
  return(status);
}



