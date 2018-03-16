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

#include "sixt.h"
#include "event.h"
#include "eventfile.h"

#define TOOLSUB ero_fits2tm_main
#include "headas_main.c"


// Number of bytes per TLM byte frame
#define N_BYTES (128)


struct Parameters {
  char EventList[MAXFILENAME];
  char OutputFile[MAXFILENAME];

  double TIMEZERO;
};


// Data structure for the byte output writing routine.
struct Binary_Output {
  unsigned char *bytes; // Byte output buffer
  FILE *fptr;           // Pointer to output file
  int n_bytes;          // Counting the bytes already written to the buffer
  int max_bytes;        // Maximum number of bytes in one TLM record
  unsigned char framecounter; // Counter for the records (NOT detector frames)
  unsigned char framenumber;  // Counter for the records belonging to one and
                              // the same detector frame
};


// Clear the given byte buffer with 'length' bytes.
void binary_output_clear_bytes(unsigned char *bytes, const int length) 
{
  int count;
  for (count=0; count<length; count++) {
    bytes[count]=0;
  }
}


// Constructor of Binary_Output
struct Binary_Output *get_Binary_Output(const int max_bytes, FILE* const fptr) 
{
  // Get memory for the object itself:
  struct Binary_Output *binary_output=
    (struct Binary_Output *)malloc(sizeof(struct Binary_Output));
  if (binary_output==NULL) return(NULL);

  // Get memory for the byte output buffer and clear it:
  binary_output->bytes=(unsigned char *)malloc(max_bytes*sizeof(unsigned char));
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
void free_Binary_Output(struct Binary_Output *binary_output) 
{
  // Free all allocated memory:
  if (NULL!=binary_output) {
    if (binary_output->bytes) {
      free(binary_output->bytes);
    }
    free(binary_output);
  }
  binary_output=NULL;
}


// This routine finishes a record and writes all bytes to the output file.
int binary_output_erosita_finish_record(struct Binary_Output* const binary_output)
{
  // If the record already contains some elements but is not full,
  // fill it with 0x00 entries.
  if (binary_output->n_bytes>12) {
    int count;
    for (count=binary_output->n_bytes; count<binary_output->max_bytes; count++) {
      binary_output->bytes[count]=0x00;
    }

    // Write the header of the record:
    // syncwords
    binary_output->bytes[0]=0xEB;
    binary_output->bytes[1]=0x90;
    // frametime MSB
    binary_output->bytes[2]=0x00;  // TODO
    binary_output->bytes[3]=0x00;
    binary_output->bytes[4]=0x00;
    binary_output->bytes[5]=0x00;
    // framecounter
    binary_output->bytes[6]=binary_output->framecounter++;
    // header
    binary_output->bytes[7]=0x1F;
    // DSP frame ID
    binary_output->bytes[8]=0x20;
    // # of 16bit elements
    binary_output->bytes[9]=(unsigned char)((binary_output->n_bytes-12)/2);
    // # of total frame   ???????????  // TODO
    binary_output->bytes[10]=0x00;
    // actual frame #
    binary_output->bytes[11]=binary_output->framenumber;

    // Write the byte buffer to the output file:
    int nbytes=fwrite(binary_output->bytes, 1, 
		      binary_output->max_bytes, binary_output->fptr);
    if (nbytes!=binary_output->max_bytes) return(-1);
    
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
				       const Event* const event)
{
  // Check for overflows:
  assert(event->pha<=0x3FFF);
  assert(event->rawx<=0xFF);
  assert(event->rawy<=0xFF);
  
  // Write the data of the event to the byte output buffer:
  binary_output->bytes[binary_output->n_bytes++]=
    0x3F & (unsigned char)(event->pha>>8);
  binary_output->bytes[binary_output->n_bytes++]=
    0xFF & (unsigned char)event->pha;
  binary_output->bytes[binary_output->n_bytes++]=
    0xFF & (unsigned char)event->rawy;
  binary_output->bytes[binary_output->n_bytes++]=
    0xFF & (unsigned char)event->rawx;

  // Check if record is full:
  if (binary_output->n_bytes>=binary_output->max_bytes) {
    if (binary_output_erosita_finish_record(binary_output)) return(-1);
  }

  return(0);
}


// This routine finishes a detector frame by adding the 3 last 4-byte elements
// containing the frame time, the number of discarded pixels and a spare element.
int binary_output_erosita_finish_frame(struct Binary_Output *binary_output, 
				       const double time) 
{
  // Write the TIME element (readout time of last detector frame):
  long ltime=(long)(time*40);
  if (ltime > 0x3FFFFFFF) {
    SIXT_ERROR("time exceeds maximum allowed value");
    return(-1);
  }
  binary_output->bytes[binary_output->n_bytes++]=
    0x40 + (unsigned char)(ltime>>24);
  binary_output->bytes[binary_output->n_bytes++]=
    (unsigned char)(ltime>>16);
  binary_output->bytes[binary_output->n_bytes++]=
    (unsigned char)(ltime>>8);
  binary_output->bytes[binary_output->n_bytes++]=
    (unsigned char)ltime;

  // Check if record is full:
  if (binary_output->n_bytes>=binary_output->max_bytes) {
    if (binary_output_erosita_finish_record(binary_output)) return(-1);
  }
  
  // Write number of DISCARDED PIXELS (assumed to be 0):
  binary_output->bytes[binary_output->n_bytes++]=0xBE;
  binary_output->bytes[binary_output->n_bytes++]=0x00;
  binary_output->bytes[binary_output->n_bytes++]=0x00;
  binary_output->bytes[binary_output->n_bytes++]=0x00;

  // Check if record is full:
  if (binary_output->n_bytes>=binary_output->max_bytes) {
    if (binary_output_erosita_finish_record(binary_output)) return(-1);
  }
  
  // Write SPARE element:
  binary_output->bytes[binary_output->n_bytes++]=0xC0;
  binary_output->bytes[binary_output->n_bytes++]=0x00;
  binary_output->bytes[binary_output->n_bytes++]=0x00;
  binary_output->bytes[binary_output->n_bytes++]=0x00;

  // Finish record in any case, even if it is not full yet:
  if (binary_output_erosita_finish_record(binary_output)) {
    return(-1);
  } else {
    binary_output->framenumber=0;
    return(0);
  }
}


int ero_fits2tm_main()
{
  struct Parameters par;

  /** FITS file containing the input event list. */
  EventFile* elf=NULL; 
  /** Output file. */
  FILE *of=NULL;

  /** Output buffer. */
  Event* eventlist=NULL;
  struct Binary_Output *binary_output=NULL;

  int status=EXIT_SUCCESS;


  // HEATOOLs: register program
  set_toolname("ero_fits2tm");
  set_toolversion("0.05");


  do { // Beginning of ERROR handling loop

    // --- Initialization ---
    char* sbuffer=NULL;

    status=ape_trad_query_file_name("EventList", &sbuffer);
    if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("failed reading the name of the input event list");
      break;
    }
    strcpy(par.EventList, sbuffer);
    free(sbuffer);

    status=ape_trad_query_file_name("OutputFile", &sbuffer);
    if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("failed reading the name of the output file");
      break;
    }
    strcpy(par.OutputFile, sbuffer);
    free(sbuffer);

    status=ape_trad_query_double("TIMEZERO", &par.TIMEZERO);
    if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("failed reading TIMEZERO");
      break;
    }

    // Open the event list FITS file.
    elf=openEventFile(par.EventList, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Check if the input file contains single-pixel events.
    char evtype[MAXMSG], comment[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "EVTYPE", evtype, comment, &status);
    if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("could not read FITS keyword 'EVTYPE'");
      break;
    }
    strtoupper(evtype);
    if (0!=strcmp(evtype, "PIXEL")) {
      status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "event type of input file is '%s' (must be 'PIXEL')", evtype);
      SIXT_ERROR(msg);
      break;
    }

    // Open the binary file for output.
    of=fopen(par.OutputFile, "w+");
    if (of==NULL) {
      status=EXIT_FAILURE;
      char msg[MAXMSG]; // Error message buffer.
      sprintf(msg, "output file '%s' could not be opened", par.OutputFile);
      SIXT_ERROR(msg);
      break;
    }

    binary_output=get_Binary_Output(N_BYTES, of);
    eventlist=(Event*)malloc(10000*sizeof(Event));
    if ((NULL==binary_output)||(NULL==eventlist)) {
      status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation failed");
      break;
    }
    // Clear the event buffer:
    int count;
    for (count=0; count<10000; count++) {
      eventlist[count].time=0.;
      eventlist[count].pha=0;
      eventlist[count].rawx=0;
      eventlist[count].rawy=0;
      eventlist[count].frame=0;
    }

    // --- END of Initialization ---


    // --- Beginning of EVENT PROCESSING ---

    // Loop over all events in the FITS file:
    headas_chat(3, "processing events ...\n");
    
    // Loop over all entries in the event list:
    int n_buffered_events=0;
    long row;
    for (row=0; row<elf->nrows; row++) {

      // Read the event from the FITS file.
      getEventFromFile(elf, row+1, &(eventlist[n_buffered_events]), &status);
      CHECK_STATUS_BREAK(status);

      // Subtract the time offset.
      eventlist[n_buffered_events].time-=par.TIMEZERO;

      if (eventlist[n_buffered_events].frame>eventlist[0].frame) {
	// Write the events to the binary output.
	for (count=0; count<n_buffered_events; count++) {
	  if (binary_output_erosita_insert_event(binary_output, &(eventlist[count]))) {
	    status=EXIT_FAILURE;
	    SIXT_ERROR("generation of binary format failed");
	    break;
	  }
	}
	// Finish the TLM record, even if it is not full yet.
	if (binary_output_erosita_finish_frame(binary_output, eventlist[0].time)) {
	  status=EXIT_FAILURE;
	  SIXT_ERROR("generation of binary format failed");
	  break;
	}
	  
	// New buffering period has started.
	eventlist[0]=eventlist[n_buffered_events];
	n_buffered_events=0;

      } // END of loop over all buffered events
      CHECK_STATUS_BREAK(status);

      n_buffered_events++;
	
    } // END of loop over all entries in the event list.
    CHECK_STATUS_BREAK(status);

    // Write the events to the binary output.
    for (count=0; count<n_buffered_events; count++) {
      if (binary_output_erosita_insert_event(binary_output, &(eventlist[count]))) {
	status=EXIT_FAILURE;
	SIXT_ERROR("generation of binary format failed");
	break;
      }
    }
    CHECK_STATUS_BREAK(status);

    if (binary_output_erosita_finish_frame(binary_output, eventlist[0].time)) {
      status=EXIT_FAILURE;
      SIXT_ERROR("generation of binary format failed");
      break;
    }
  
  } while (0); // END of ERROR handling loop


  // --- Clean up ---

  // Close files
  if (NULL!=of) fclose(of);
  freeEventFile(&elf, &status);

  // Release memory.
  free_Binary_Output(binary_output);
  if (NULL!=eventlist) free(eventlist);
  
  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}

