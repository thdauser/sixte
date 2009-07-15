#include "sixt.h"
#include "erositaevent.h"


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


