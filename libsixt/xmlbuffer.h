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

/** Data structure for include handling */
struct XMLIncludeHandler {

  /** Flag if the preprocessed XMLBuffer contained any further includes
      to be expanded. */
  int further_includes;

  /** Output buffer for included XML data. */
  struct XMLBuffer* include_buffer;

  /** Output buffer for processed XML data. */
  struct XMLBuffer* output_buffer;

  /** filepath of the xml file */
  char xmlfile[MAXFILENAME];

  int status;
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

/** Expand the included XML files in the GenDet XML
    description. */
void expandIncludesXML(struct XMLBuffer* const buffer, 
		       const char* filename, 
		       int* const status);

/** Expand the loops and arithmetic operations in the GenDet XML
    description. */
void expandXML(struct XMLBuffer* const buffer, int* const status);


/** Read the string value of an XML element. */
void getXMLAttributeString(const char** attr, const char* const key, 
			   char* const value);

/** Read the float value of an XML element. */
float getXMLAttributeFloat(const char** attr, const char* const key);

/** Read the double value of an XML element. */
double getXMLAttributeDouble(const char** attr, const char* const key);

/** Read the integer value of an XML element. */
int getXMLAttributeInt(const char** attr, const char* const key);

/** Read the long value of an XML element. */
long getXMLAttributeLong(const char** attr, const char* const key);


#endif /* XMLBUFFER_H */
