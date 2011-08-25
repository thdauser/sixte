#ifndef XMLBUFFER_H 
#define XMLBUFFER_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Buffer for XML code read from the file and expanded in order to
    handle loops. */
struct XMLBuffer {
  char* text;
  unsigned long maxlength;
};

/** Data structure given to the XML Pre-Parser. */
struct XMLPreParseData {

  /** Flag if the preprocessed XMLBuffer contained any further loops
      to be expanded. */
  int further_loops;

  /** Current loop depth. */
  int loop_depth;
  /** Start, end, and increment of the outermost loop. */
  int loop_start, loop_end, loop_increment;
  /** Loop counter variable. This variable can be used in the XML text
      as $[NAME]. */
  char loop_variable[MAXMSG];
  
  /** Output buffer for processed XML data. */
  struct XMLBuffer* output_buffer;
  /** Buffer for XML code inside the loop. */
  struct XMLBuffer* loop_buffer;
  
  int status;
};


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor of XMLBuffer. */
struct XMLBuffer* newXMLBuffer(int* const status);

/** Destructor of XMLBuffer. Release the memory from the string
    buffer. */
void freeXMLBuffer(struct XMLBuffer** const buffer);

/** Add a string to the XMLBuffer. If the buffer size is to small,
    allocate additional memory. */
void addString2XMLBuffer(struct XMLBuffer* const buffer, 
			 const char* const string,
			 int* const status);

/** Expand the loops and arithmetic operations in the GenDet XML
    description. */
void expandXML(struct XMLBuffer* const buffer, int* const status);


#endif /* XMLBUFFER_H */
