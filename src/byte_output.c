#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <ctype.h>
#include <math.h>

#include "eventlist.types.h"

// Data structure for the byte output writing routine.
struct Byte_Output {
  unsigned char *bytes;  // Byte output buffer
  FILE *fptr;            // Pointer to output file
  int n_bytes;           // Counting the bytes already written to the buffer
  int max_bytes;         // Maximum number of bytes in one TLM record
  unsigned char framecounter;  // Counter for the records (NOT detector frames)
  unsigned char framenumber;   // Counter for the records belonging to one and
                               // the same detector frame
};




// Clear the given byte buffer with 'length' bytes.
void byte_output_clear_bytes(unsigned char *bytes, const int length) {
  int count;
  for (count=0; count<length; count++) {
    bytes[count] = 0;
  }
}


// Constructor of Byte_Output
struct Byte_Output *get_Byte_Output(const int max_bytes, FILE *fptr) {
  // Get memory for the object itself:
  struct Byte_Output *byte_output = 
    (struct Byte_Output *) malloc(sizeof(struct Byte_Output));
  if (byte_output==NULL) return(NULL);

  // Get memory for the byte output buffer and clear it:
  byte_output->bytes = (unsigned char *) malloc(max_bytes*sizeof(unsigned char));
  if (byte_output->bytes==NULL) return(NULL);
  byte_output_clear_bytes(byte_output->bytes, max_bytes);

  // Set object configuration:
  byte_output->n_bytes=12; // -> byte address of the first element in the record
  byte_output->max_bytes=max_bytes;
  byte_output->framecounter=0;
  byte_output->framenumber=0;

  return(byte_output);
}


// Destructor of Byte_Output
void free_Byte_Output(struct Byte_Output *byte_output) {
  // Free all allocated memory:
  if (byte_output) {
    if (byte_output->bytes) {
      free(byte_output->bytes);
    }
    free(byte_output);
  }
  byte_output=NULL;
}


// Routine which is called to write an eROSITA event to the Byte_Output.
// Return value is '0' if everything is ok, otherwise the function returns '-1'.
int byte_output_erosita_insert_event(struct Byte_Output *byte_output,
				     struct Event *event)
{
  // Check for overflows:
  if (event->pha > 0x3FFF) return (-1);

  // Write the data of the event to the byte output buffer:
  byte_output->bytes[byte_output->n_bytes++] = 
    0x3F & (unsigned char)(event->pha>>8);
  byte_output->bytes[byte_output->n_bytes++] = 
    0xFF & (unsigned char)event->pha;
  byte_output->bytes[byte_output->n_bytes++] = 
    0xFF & (unsigned char)event->yi;
  byte_output->bytes[byte_output->n_bytes++] = 
    0xFF & (unsigned char)event->xi;

  // Check if record is full:
  if (byte_output->n_bytes>byte_output->max_bytes) {
    if (byte_output_erosita_finish_record(byte_output)) return(-1);
  }

  return(0);
}


// This routine finihes a record and writes all bytes to the output file.
int byte_output_erosita_finish_record(struct Byte_Output *byte_output){
  int count;

  // If the record already contains some elements but is not full,
  // fill it with 0x00 entries.
  if (byte_output->n_bytes>12) {
    for (count=byte_output->n_bytes; count<=byte_output->max_bytes; count++) {
      byte_output->bytes[count] = 0x00;
    }

    // Write the header of the record:
    // syncwords
    byte_output->bytes[0] = 0xEB;
    byte_output->bytes[1] = 0x90;
    // frametime MSB
    byte_output->bytes[2] = 0x00;
    byte_output->bytes[3] = 0x00;
    byte_output->bytes[4] = 0x00;
    byte_output->bytes[5] = 0x00;
    // framecounter
    byte_output->bytes[6] = byte_output->framecounter++;
    // header
    byte_output->bytes[7] = 0x1F;
    // DSP frame ID
    byte_output->bytes[8] = 0x20;
    // # of 16bit elements
    byte_output->bytes[9] = (unsigned char)((byte_output->n_bytes-12)/2);
    // # of total frame   ???????????
    byte_output->bytes[10] = 0x00;
    // actual frame #
    byte_output->bytes[11] = byte_output->framenumber;


    // Write the byte buffer to the output file:
    int nbytes = fwrite (byte_output->bytes, 1, 
			 byte_output->max_bytes, byte_output->fptr);
    if (nbytes != byte_output->max_bytes) return(-1);
    

    // Reset variables:
    byte_output->n_bytes=12;

    // Clear output buffer:
    byte_output_clear_bytes(byte_output->bytes, byte_output->max_bytes);
  } // END (byte_output->n_bytes>12)

  return(0);
}



// This routine finihes a detector frame by adding the 3 last 4-byte elements
// containing the frame time, the number of discarded pixels and a spare element.
int byte_output_erosita_finish_frame(struct Byte_Output *byte_output, 
				     const double time) 
{
  // Write the TIME element (readout time of last detector frame):
  long ltime = (long)(time*40);
  if (ltime > 0x3FFFFFFF)  return(-1);
  byte_output->bytes[byte_output->n_bytes++] = 
    0x40 + (unsigned char)(ltime>>24);
  byte_output->bytes[byte_output->n_bytes++] = 
    (unsigned char)(ltime>>16);
  byte_output->bytes[byte_output->n_bytes++] = 
    (unsigned char)(ltime>>8);
  byte_output->bytes[byte_output->n_bytes++] = 
    (unsigned char)ltime;

  // Check if record is full:
  if (byte_output->n_bytes>byte_output->max_bytes) {
    if (byte_output_erosita_finish_record(byte_output)) return(-1);
  }

  
  // Write number of DISCARDED PIXELS (assumed to be 0):
  byte_output->bytes[byte_output->n_bytes++] = 0xBE;
  byte_output->bytes[byte_output->n_bytes++] = 0x00;
  byte_output->bytes[byte_output->n_bytes++] = 0x00;
  byte_output->bytes[byte_output->n_bytes++] = 0x00;

  // Check if record is full:
  if (byte_output->n_bytes>byte_output->max_bytes) {
    if (byte_output_erosita_finish_record(byte_output)) return(-1);
  }
  

  // Write SPARE element:
  byte_output->bytes[byte_output->n_bytes++] = 0xC0;
  byte_output->bytes[byte_output->n_bytes++] = 0x00;
  byte_output->bytes[byte_output->n_bytes++] = 0x00;
  byte_output->bytes[byte_output->n_bytes++] = 0x00;

  // Finish record in any case, even if it is not full yet:
  if (byte_output_erosita_finish_record(byte_output)) {
    return(-1);
  } else {
    byte_output->framenumber=0;
    return(0);
  }
}


