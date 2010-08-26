#include "gendet.h"

////////////////////////////////////////////////////////////////////
// Static function declarations
////////////////////////////////////////////////////////////////////

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


/** Return the index of the detector line affected by the specified
    y-value. If the y-value is outside the line region, the return
    value is -1. */
static int getGenDetAffectedLine(GenDet* const det, const double y);

/** Return the index of the detector column affected by the specified
    x-value. If the x-value is outside the line region, the return
    value is -1. */
static int getGenDetAffectedColumn(GenDet* const det, const double x);

/** Return the index of the bin affected by the specified
    x-position. The bin grid is defined by the following WCS compliant
    values: reference pixel (rpix), reference value (rval), and pixel
    delta (delt). If the specified x-value is outside the bins, the
    function return value is -1. */
static int getAffectedIndex(const double x, const float rpix, 
			    const float rval, const float delt, 
			    const int width);


////////////////////////////////////////////////////////////////////
// Code
////////////////////////////////////////////////////////////////////

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
  det->xwidth=-1;
  det->ywidth=-1;
  det->xrpix =-1.;
  det->yrpix =-1.;
  det->xrval =-1.;
  det->yrval =-1.;
  det->xdelt =-1.;
  det->ydelt =-1.;
  strcpy(det->rmf_filename, "");
  det->readout_trigger = 0;

  // Read in the XML definition of the detector.
  parseGenDetXML(det, filename, status);
  if (EXIT_SUCCESS!=*status) return(det);

  // Check if all required parameters have been read successfully from 
  // the XML file.
  if (-1==det->xwidth) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for x-width of GenDet pixel array!\n", 
		   *status);
    return(det);    
  }  
  if (-1==det->ywidth) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for y-width of GenDet pixel array!\n", 
		   *status);
    return(det);    
  }

  if (0>det->xrpix) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for x reference pixel of GenDet!\n", 
		   *status);
    return(det);    
  }
  if (0>det->yrpix) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for y reference pixel of GenDet!\n", 
		   *status);
    return(det);    
  }

  if (0>det->xrval) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for x reference value of GenDet!\n", 
		   *status);
    return(det);    
  }
  if (0>det->yrval) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for y reference value of GenDet!\n", 
		   *status);
    return(det);    
  }

  if (0>det->xdelt) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for x pixel width of GenDet!\n", 
		   *status);
    return(det);    
  }
  if (0>det->ydelt) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for y pixel width of GenDet!\n", 
		   *status);
    return(det);    
  }
  
  if (0==strlen(det->rmf_filename)) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for response file of GenDet!\n", 
		   *status);
    return(det);    
  }

  if (0==det->readout_trigger) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for the readout trigger of GenDet!\n", 
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

  // Load the detector response file.
  det->rmf = loadRMF(det->rmf_filename, status);
  if (EXIT_SUCCESS!=*status) return(det);

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
  headas_chat(5, "read detector setup from XML file '%s' ...\n", filename);

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
  char Uvalue[MAXMSG];     // Upper case version of XML attribute value

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

    else if (!strcmp(Uelement, "WCS")) {
      if (!strcmp(Uattribute, "XRPIX")) {
	xmldata->det->xrpix = (float)atof(attr[i+1]);
      } else if (!strcmp(Uattribute, "YRPIX")) {
	xmldata->det->yrpix = (float)atof(attr[i+1]);
      } else if (!strcmp(Uattribute, "XRVAL")) {
	xmldata->det->xrval = (float)atof(attr[i+1]);
      } else if (!strcmp(Uattribute, "YRVAL")) {
	xmldata->det->yrval = (float)atof(attr[i+1]);
      } else if (!strcmp(Uattribute, "XDELT")) {
	xmldata->det->xdelt = (float)atof(attr[i+1]);
      } else if (!strcmp(Uattribute, "YDELT")) {
	xmldata->det->ydelt = (float)atof(attr[i+1]);
      }
    }
    
    else if (!strcmp(Uelement, "RESPONSE")) {
      if (!strcmp(Uattribute, "FILENAME")) {
	strcpy(xmldata->det->rmf_filename, attr[i+1]);
      }
    }

    else if (!strcmp(Uelement, "READOUT")) {
      if (!strcmp(Uattribute, "MODE")) {
	strcpy(Uvalue, attr[i+1]);
	strtoupper(Uvalue);
	if (!strcmp(Uvalue, "TIME")) {
	  xmldata->det->readout_trigger = GENDET_TIME_TRIGGERED;
	} else if (!strcmp(Uvalue, "EVENT")) {
	  xmldata->det->readout_trigger = GENDET_EVENT_TRIGGERED;
	}
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



static int getGenDetAffectedLine(GenDet* const det, const double y)
{
  return(getAffectedIndex(y, det->yrpix, det->yrval, det->ydelt, det->ywidth));
}



static int getGenDetAffectedColumn(GenDet* const det, const double x)
{
  return(getAffectedIndex(x, det->xrpix, det->xrval, det->xdelt, det->xwidth));
}



static int getAffectedIndex(const double x, const float rpix, 
			    const float rval, const float delt, 
			    const int width)
{
  int index = ((int)((x+ (rpix-0.5)*delt -rval)/delt +1.))-1;
  //                  avoid (int)(-0.5) = 0     <-----|----|
  if (index>=width) index=-1;
  return(index);
}



void addGenDetPhotonImpact(GenDet* const det, const Impact* const impact, int* const status)
{
  // Determine the affected detector line and column.
  int line   = getGenDetAffectedLine  (det, impact->position.y);
  int column = getGenDetAffectedColumn(det, impact->position.x);

  printf("%e %d\n", impact->position.x, column);
  printf("%e %d\n", impact->position.y, line);

  // Check if the returned values are valid line and column indices.
  if ((0>column) || (0>line)) {
    return;
  }

}

