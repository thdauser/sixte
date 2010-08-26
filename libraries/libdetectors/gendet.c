#include "gendet.h"


/** Parse the GenDet definition from an XML file. */
static void parseGenDetXML(GenDet* const det, const char* const filename, int* const status);

/** Handler for the start of an XML element. */
static void GenDetXMLElementStart(void *data, const char *el, const char **attr);
/** Handler for the end of an XML element. */
static void GenDetXMLElementEnd(void *data, const char *el);

/** Data structure given to the XML handler to transfer data. */
struct XMLData {
  GenDet* det;
  int depth;
};


GenDet* newGenDet(const char* const filename, int* const status) {
  GenDet* det=NULL;

  // Allocate memory.
  det=(GenDet*)malloc(sizeof(GenDet));
  if (NULL==det) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenDet failed!\n", *status);
    return(det);
  }

  // Initialize all pointers with NULL.
  det->line=NULL;
  det->rmf =NULL;

  // Set initial values before parsing the parameters from the XML file.
  det->xwidth=0;
  det->ywidth=0;

  // Read in the XML definition of the detector.
  parseGenDetXML(det, filename, status);
  if (EXIT_SUCCESS!=*status) return(det);

  // Check if all required parameters have been read successfully from 
  // the XML file.
  if (0==det->xwidth) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for x-width of GenDet pixel array!\n", 
		   *status);
    return(det);    
  }  
  if (0==det->ywidth) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for y-width of GenDet pixel array!\n", 
		   *status);
    return(det);    
  }
  // END of checking, if all detector parameters have successfully been 
  // read from the XML file.

  // Allocate memory for the pixels.
  det->line=(GenDetLine**)malloc(det->ywidth*sizeof(GenDetLine*));
  if (NULL==det->line) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenDet failed!\n", *status);
    return(det);
  }
  int i;
  for (i=0; i<det->ywidth; i++) {
    det->line[i] = newGenDetLine(det->xwidth, status);
    if (EXIT_SUCCESS!=*status) return(det);
  }

  return(det);
}



void destroyGenDet(GenDet** det)
{
  if (NULL!=*det) {
    if (NULL!=(*det)->line) {
      int i;
      for (i=0; i<(*det)->ywidth; i++) {
	destroyGenDetLine(&(*det)->line[i]);
      }
      free((*det)->line);
    }

    free(*det);
    *det=NULL;
  }
}



static void parseGenDetXML(GenDet* const det, const char* const filename, int* const status)
{
  // Open the specified file.
  FILE* xmlfile = fopen(filename, "r");
  if (NULL==xmlfile) {
    *status = EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "Error: Failed opening GenDet definition XML "
	    "file '%s' for read access!\n", filename);
    HD_ERROR_THROW(msg, *status);
    return;
  }

  // Parse the XML data from the file using the expat library.

  // Get an XML_Parser object.
  XML_Parser parser = XML_ParserCreate(NULL);
  if (NULL==parser) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for XML parser!\n", *status);
    return;
  }

  // Set data that is passed to the handler functions.
  struct XMLData xmldata = {
    .depth = 0,
    .det   = det
  };
  XML_SetUserData(parser, &xmldata);

  // Set the handler functions.
  XML_SetElementHandler(parser, GenDetXMLElementStart, GenDetXMLElementEnd);

  int done =0;
  int len; // Number of chars in buffer.
  const int buffer_size=10;
  // Input buffer with an additional byte at the end for the 
  // termination of the string.
  char buffer[buffer_size+1];

  do {
    // Get a piece of input into the buffer.
    len = fread(buffer, 1, buffer_size, xmlfile);
    buffer[len]='\0'; // Terminate the string.
    done = feof(xmlfile);
    //    printf("%s", buffer);//RM
    if (!XML_Parse(parser, buffer, len, done)) {
      // Parse error.
      *status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "Error: Parsing XML file '%s' failed at line %d:\n%s\n", 
	      filename,
	      (int)XML_GetCurrentLineNumber(parser),
	      XML_ErrorString(XML_GetErrorCode(parser)));
      HD_ERROR_THROW(msg, *status);
      return;
    }
  } while (!done);
  XML_ParserFree(parser);

  // Close the file handler.
  fclose(xmlfile);
}



static void GenDetXMLElementStart(void *data, const char *el, const char **attr) 
{
  struct XMLData* xmldata = (struct XMLData*)data;
  char Uelement[MAXMSG];   // Upper case version of XML element
  char Uattribute[MAXMSG]; // Upper case version of XML attribute

  // Convert the element to an upper case string.
  strcpy(Uelement, el);
  strtoupper(Uelement);

  // Loop over the different attributes.
  int i;
  for (i=0; attr[i]; i+=2) {

    // Convert the attribute to an upper case string.
    strcpy(Uattribute, attr[i]);
    strtoupper(Uattribute);

    // Check the XML element name.
    if (!strcmp(Uelement, "DIMENSIONS")) {
      if (!strcmp(Uattribute, "XWIDTH")) {
	xmldata->det->xwidth = atoi(attr[i+1]);
      } else if (!strcmp(Uattribute, "YWIDTH")) {
	xmldata->det->ywidth = atoi(attr[i+1]);
      }
    }
    
  } 
  // END of loop over different attributes.
}



static void GenDetXMLElementEnd(void *data, const char *el) 
{
  struct XMLData* xmldata = (struct XMLData*)data;

  xmldata->depth--;

  return;
}
