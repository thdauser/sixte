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

/** Data structure given to the XML Hexagon-Parser. */
struct XMLHexParseData {
	/** Radius [m] of the hexagon to realize */
	double radius;
	/** Pitch [m] in the two directions between the pixels */
	double pixelpitch;
	/** Bool to know whether we are still in the loop */
	char inside_loop;

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

/** Expand the hexagonal detector loop in the advanced detector definition */
void expandHexagon(struct XMLBuffer* const buffer, int* const status);


#endif /* XMLBUFFER_H */
