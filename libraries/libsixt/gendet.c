#include "gendet.h"

////////////////////////////////////////////////////////////////////
// Static data type declarations
////////////////////////////////////////////////////////////////////


/** Data structure given to the XML handler to transfer data. */
struct XMLParseData {
  GenDet* det;
  int status;
};

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


////////////////////////////////////////////////////////////////////
// Static function declarations
////////////////////////////////////////////////////////////////////


/** Parse the GenDet definition from an XML file. */
static void parseGenDetXML(GenDet* const det, const char* const filename, 
			   int* const status);

/** Handler for the start of an XML element. */
static void GenDetXMLElementStart(void* data, const char* el, 
				  const char** attr);
/** Handler for the end of an XML element. */
static void GenDetXMLElementEnd(void* data, const char* el);


/** Expand the loops and arithmetic operations in the GenDet XML
    description. */
static void expandXML(struct XMLBuffer* const buffer, int* const status);
/** Handler for the start of an XML element. */
static void expandXMLElementStart(void* data, const char* el, 
				  const char** attr);
/** Handler for the end of an XML element. */
static void expandXMLElementEnd(void* data, const char* el);


/** Constructor of XMLBuffer. */
static struct XMLBuffer* newXMLBuffer(int* const status);

/** Destructor of XMLBuffer. Release the memory from the string
    buffer. */
static void destroyXMLBuffer(struct XMLBuffer** const buffer);

/** Add a string to the XMLBuffer. If the buffer size is to small,
    allocate additional memory. */
static void addString2XMLBuffer(struct XMLBuffer* const buffer, 
				const char* const string,
				int* const status);

/** Copy an XMLBuffer string from the source to the destination. */
static void copyXMLBuffer(struct XMLBuffer* const destination,
			  struct XMLBuffer* const source,
			  int* const status);

/** Replace all occurences of the string 'old' in the XMLBuffer text
    by the string 'new'. */
static void replaceInXMLBuffer(struct XMLBuffer* const buffer, 
			       const char* const old,
			       const char* const new, 
			       int* const status);

/** Execute and replace arithmetic operations in the XMLBuffer
    text. */
static void execArithmeticOpsInXMLBuffer(struct XMLBuffer* const buffer,
					 int* const status);


////////////////////////////////////////////////////////////////////
// Program Code
////////////////////////////////////////////////////////////////////


GenDet* newGenDet(const char* const filename, int* const status) 
{
  // Allocate memory.
  GenDet* det=(GenDet*)malloc(sizeof(GenDet));
  if (NULL==det) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenDet failed!\n", *status);
    return(det);
  }

  // Initialize all pointers with NULL.
  det->pixgrid=NULL;
  det->split  =NULL;
  det->line=NULL;
  det->psf =NULL;
  det->vignetting=NULL;
  det->coded_mask=NULL;
  det->arf =NULL;
  det->rmf =NULL;
  det->clocklist =NULL;
  det->badpixmap =NULL;
  det->eventfile =NULL;
  det->grading=NULL;

  // Get empty GenPixGrid.
  det->pixgrid = newGenPixGrid(status);
  if (EXIT_SUCCESS!=*status) return(det);

  // Get empty ClockList.
  det->clocklist = newClockList(status);
  if (EXIT_SUCCESS!=*status) return(det);

  // Get empty split model.
  det->split = newGenSplit(status);
  if (EXIT_SUCCESS!=*status) return(det);

  // Read in the XML definition of the detector.
  parseGenDetXML(det, filename, status);
  if (EXIT_SUCCESS!=*status) return(det);
    
  // Allocate memory for the pixels.
  det->line=(GenDetLine**)malloc(det->pixgrid->ywidth*sizeof(GenDetLine*));
  if (NULL==det->line) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenDet failed!\n", *status);
    return(det);
  }
  int ii;
  for (ii=0; ii<det->pixgrid->ywidth; ii++) {
    det->line[ii] = newGenDetLine(det->pixgrid->xwidth, status);
    if (EXIT_SUCCESS!=*status) return(det);
  }

  return(det);
}



void destroyGenDet(GenDet** const det, int* const status)
{
  if (NULL!=*det) {
    // Destroy the pixel array.
    if (NULL!=(*det)->line) {
      int i;
      for (i=0; i<(*det)->pixgrid->ywidth; i++) {
	destroyGenDetLine(&(*det)->line[i]);
      }
      free((*det)->line);
    }

    // Destroy the ClockList.
    destroyClockList(&(*det)->clocklist);

    // Destroy the GenPixGrid.
    destroyGenPixGrid(&(*det)->pixgrid);

    // Destroy the split model.
    destroyGenSplit(&(*det)->split);

    // Destroy the BadPixMap.
    destroyBadPixMap(&(*det)->badpixmap);

    // Destroy the pattern identifier.
    destroyGenEventGrading(&(*det)->grading);

    // Close the event file.
    destroyGenEventFile(&(*det)->eventfile, status);

    // Destroy the PSF.
    destroyPSF(&(*det)->psf);

    // Destroy the CodedMask.
    destroyCodedMask(&(*det)->coded_mask);

    // Destroy the vignetting Function.
    destroyVignetting(&(*det)->vignetting);

    free(*det);
    *det=NULL;
  }
}



static void addString2XMLBuffer(struct XMLBuffer* const buffer, 
				const char* const string,
				int* const status)
{
  // Check if a valid buffer is specified.
  if (NULL==buffer) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: NULL pointer to XMLBuffer!\n", *status);
    return;
  }
    
  // Check if the buffer is empty.
  if (NULL==buffer->text) {
    // Allocate memory for the first chunk of bytes.
    buffer->text=(char*)malloc((MAXMSG+1)*sizeof(char));
    if (NULL==buffer->text) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for XMLBuffer failed!\n", *status);
      return;
    }
    buffer->text[0]='\0';
    buffer->maxlength=MAXMSG;
  }

  // Check if the buffer contains sufficient memory to add the new string.
  if (strlen(buffer->text)+strlen(string)>=buffer->maxlength) {
    // Allocate the missing memory.
    int new_length = strlen(buffer->text) + strlen(string);
    buffer->text=(char*)realloc(buffer->text, (new_length+1)*sizeof(char));
    if (NULL==buffer->text) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for XMLBuffer failed!\n", *status);
      return;
    }
    buffer->maxlength=new_length;
  }

  // Append the new string to the existing buffer.
  strcat(buffer->text, string);
}



static void copyXMLBuffer(struct XMLBuffer* const destination,
			  struct XMLBuffer* const source,
			  int* const status)
{
  // Adapt memory size.
  destination->text = (char*)realloc(destination->text,
				     (source->maxlength+1)*sizeof(char));
  if (NULL==destination->text) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: memory allocation for XMLBuffer failed!\n", *status);
    return;
  }
  destination->maxlength=source->maxlength;

  // Copy content.
  strcpy(destination->text, source->text);
}



static struct XMLBuffer* newXMLBuffer(int* const status)
{
  struct XMLBuffer* buffer=(struct XMLBuffer*)malloc(sizeof(struct XMLBuffer));
  if (NULL==buffer) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for XMLBuffer failed!\n", *status);
    return(buffer);
  }

  buffer->text=NULL;
  buffer->maxlength=0;

  return(buffer);
}



static void destroyXMLBuffer(struct XMLBuffer** const buffer)
{
  if (NULL!=*buffer) {
    if (NULL!=(*buffer)->text) {
      free((*buffer)->text);
    }
    free(*buffer);
    *buffer=NULL;
  }
}



static void parseGenDetXML(GenDet* const det, const char* const filename, int* const status)
{
  headas_chat(5, "read detector setup from XML file '%s' ...\n", filename);

  // Set initial values before parsing the parameters from the XML file.
  det->pixgrid->xwidth =-1;
  det->pixgrid->ywidth =-1;
  det->pixgrid->xrpix  =-1.;
  det->pixgrid->yrpix  =-1.;
  det->pixgrid->xrval  =-1.;
  det->pixgrid->yrval  =-1.;
  det->pixgrid->xdelt  =-1.;
  det->pixgrid->ydelt  =-1.;
  det->pixgrid->xborder= 0.;
  det->pixgrid->yborder= 0.;
  det->readout_trigger = 0;
  det->cte             = 1.;
  det->threshold_readout_lo_PHA    = -1;
  det->threshold_readout_up_PHA    = -1;
  det->threshold_readout_lo_keV    =  0.;
  det->threshold_readout_up_keV    = -1.;
  det->threshold_event_lo_keV      =  0.;
  det->threshold_split_lo_keV      =  0.;
  det->threshold_split_lo_fraction =  0.;
  det->fov_diameter = 0.;
  det->focal_length = 0.;
  // Set string variables to empty strings.
  strcpy(det->eventfile_template, "");
  strcpy(det->patternfile_template, "");


  // Read the XML data from the file.
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

  // The data is read from the XML file and stored in xmlbuffer
  // without any modifications.
  struct XMLBuffer* xmlbuffer = newXMLBuffer(status);
  if (EXIT_SUCCESS!=*status) return;

  // Input buffer with an additional byte at the end for the 
  // termination of the string.
  const int buffer_size=256;
  char buffer[buffer_size+1];
  // Number of chars in buffer.
  int len;

  // Read all data from the file.
  do {
    // Get a piece of input into the buffer.
    len = fread(buffer, 1, buffer_size, xmlfile);
    buffer[len]='\0'; // Terminate the string.
    addString2XMLBuffer(xmlbuffer, buffer, status);
    if (EXIT_SUCCESS!=*status) return;
  } while (!feof(xmlfile));

  // Close the file handler to the XML file.
  fclose(xmlfile);



  // Before acutally parsing the XML code, expand the loops and 
  // arithmetic operations in the GenDet XML description.
  // The expansion algorithm repeatetly scans the XML code and
  // searches for loop tags. It replaces the loop tags by repeating
  // the contained XML code.
  expandXML(xmlbuffer, status);
  if (EXIT_SUCCESS!=*status) return;



  // Parse XML code in the xmlbuffer using the expat library.
  // Get an XML_Parser object.
  XML_Parser parser = XML_ParserCreate(NULL);
  if (NULL==parser) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for XML parser!\n", *status);
    return;
  }

  // Set data that is passed to the handler functions.
  struct XMLParseData xmlparsedata = {
    .det   = det,
    .status = EXIT_SUCCESS
  };
  XML_SetUserData(parser, &xmlparsedata);

  // Set the handler functions.
  XML_SetElementHandler(parser, GenDetXMLElementStart, GenDetXMLElementEnd);

  // Parse all the data in the string buffer.
  const int done=1;
  if (!XML_Parse(parser, xmlbuffer->text, strlen(xmlbuffer->text), done)) {
    // Parse error.
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "Error: Parsing XML file '%s' failed:\n%s\n", 
	    filename, XML_ErrorString(XML_GetErrorCode(parser)));
    printf("%s", xmlbuffer->text);
    HD_ERROR_THROW(msg, *status);
    return;
  }
  // Check for errors.
  if (EXIT_SUCCESS!=xmlparsedata.status) {
    *status = xmlparsedata.status;
    return;
  }
  XML_ParserFree(parser);



  // Remove the XML string buffer.
  destroyXMLBuffer(&xmlbuffer);



  // Check if all required parameters have been read successfully from 
  // the XML file.
  if (-1==det->pixgrid->xwidth) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("no specification for x-width of GenDet pixel array in XML file");
    return;
  }  
  if (-1==det->pixgrid->ywidth) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("no specification for y-width of GenDet pixel array in XML file");
    return;
  }

  if (0>det->pixgrid->xrpix) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("no specification for x reference pixel of GenDet in XML file");
    return;    
  }
  if (0>det->pixgrid->yrpix) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("no specification for y reference pixel of GenDet in XML file");
    return;    
  }

  if (0>det->pixgrid->xrval) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("no specification for x reference value of GenDet in XML file");
    return;    
  }
  if (0>det->pixgrid->yrval) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("no specification for y reference value of GenDet in XML file");
    return;    
  }

  if (0>det->pixgrid->xdelt) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("no specification for x pixel width of GenDet in XML file");
    return;    
  }
  if (0>det->pixgrid->ydelt) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("no specification for y pixel width of GenDet in XML file");
    return;    
  }
  
  if (0.>det->pixgrid->xborder) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("invalid specification for x-border of pixels in XML file");
    return;    
  }
  if (0.>det->pixgrid->yborder) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("invalid specification for y-border of pixels in XML file");
    return;    
  }

  if (NULL==det->rmf) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("no specification for response file (RMF/RSP) in XML file");
    return;    
  }

  if (NULL==det->arf) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("no specification for ARF in XML file");
    return;    
  }

  if (0.>=det->focal_length) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("no specification for the focal length of the telescope "
	       "in the XML file");
    return;    
  }
  if (0.>=det->fov_diameter) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("no specification for the diameter of the telescope "
	       "FoV in the XML file");
    return;    
  }

  if (0==det->readout_trigger) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("no specification for the readout trigger of GenDet in the XML file");
    return;
  }

  if (0==strlen(det->eventfile_template)) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("no event file template specified in the XML file");
    return;    
  }

  if (0==strlen(det->patternfile_template)) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("no pattern file template specified in the XML file");
    return;    
  }

  if (GS_NONE!=det->split->type) {
    if (det->split->par1<=0.) {
      *status = EXIT_FAILURE;
      SIXT_ERROR("no valid split model parameters in the XML file");
      return;    
    }
  }
  // END of checking, if all detector parameters have successfully been 
  // read from the XML file.

  // If any thresholds have been specified in terms of PHA value,
  // set the corresponding charge threshold to the [keV] values
  // according to the RMF. If a charge threshold is given in addition,
  // its value is overwritten by the charge corresponding to the PHA 
  // specification. I.e., the PHA thresholds have a higher priority.
  if (det->threshold_readout_lo_PHA>-1) {
    det->threshold_readout_lo_keV = 
      getEBOUNDSEnergy(det->threshold_readout_lo_PHA, det->rmf, -1);
    headas_chat(3, "set lower readout threshold to %.3lf keV (PHA %ld)\n", 
		det->threshold_readout_lo_keV, det->threshold_readout_lo_PHA);
  }
  if (det->threshold_readout_up_PHA>-1) {
    det->threshold_readout_up_keV = 
      getEBOUNDSEnergy(det->threshold_readout_up_PHA, det->rmf,  1);
    headas_chat(3, "set upper readout threshold to %.3lf keV (PHA %ld)\n", 
		det->threshold_readout_up_keV, det->threshold_readout_up_PHA);
  }
}



static void getAttribute(const char** attr, const char* const key, char* const value)
{
  char Uattribute[MAXMSG]; // Upper case version of XML attribute
  char Ukey[MAXMSG];       // Upper case version of search expression

  // Convert the search expression to an upper case string.
  strcpy(Ukey, key);
  strtoupper(Ukey);

  int i;
  for (i=0; attr[i]; i+=2) {  
    // Convert the attribute to an upper case string.
    strcpy(Uattribute, attr[i]);
    strtoupper(Uattribute);
    if (!strcmp(Uattribute, Ukey)) {
      strcpy(value, attr[i+1]);
      return;
    }
  }
  // Keyword was not found
  strcpy(value, "");
  return;
}



static void GenDetXMLElementStart(void* parsedata, const char* el, const char** attr) 
{
  struct XMLParseData* xmlparsedata = (struct XMLParseData*)parsedata;
  char Uelement[MAXMSG];   // Upper case version of XML element
  char Uattribute[MAXMSG]; // Upper case version of XML attribute
  char Uvalue[MAXMSG];     // Upper case version of XML attribute value

  // Check if an error has occurred previously.
  if (EXIT_SUCCESS!=xmlparsedata->status) return;

  // Convert the element to an upper case string.
  strcpy(Uelement, el);
  strtoupper(Uelement);

  // Elements without attributes.
  if (!strcmp(Uelement, "LINESHIFT")) {
    CLLineShift* cllineshift = newCLLineShift(&xmlparsedata->status);
    append2ClockList(xmlparsedata->det->clocklist, CL_LINESHIFT, 
		     cllineshift, &xmlparsedata->status);
    
  } else { 
    
    // Elements with attributes.

    if (!strcmp(Uelement, "READOUTLINE")) {
      char buffer[MAXMSG]; // String buffer.
      getAttribute(attr, "LINEINDEX", buffer);
      int lineindex    = atoi(buffer);
      if (lineindex<0) {
	xmlparsedata->status=EXIT_FAILURE;
	HD_ERROR_THROW("Error: Negative index for readout line!\n", xmlparsedata->status);
	return;
      }
      getAttribute(attr, "READOUTINDEX", buffer);
      int readoutindex = atoi(buffer);
      if (readoutindex<0) {
	xmlparsedata->status=EXIT_FAILURE;
	HD_ERROR_THROW("Error: Negative index for readout line!\n", xmlparsedata->status);
	return;
      }
      CLReadoutLine* clreadoutline = newCLReadoutLine(lineindex,
						      readoutindex,
						      &xmlparsedata->status);
      append2ClockList(xmlparsedata->det->clocklist, CL_READOUTLINE, 
		       clreadoutline, &xmlparsedata->status);
	
    } else if (!strcmp(Uelement, "EVENTGRADING")) {
      char buffer[MAXMSG]; // String buffer.
      getAttribute(attr, "INVALID", buffer);
      int invalid       = atoi(buffer);
      getAttribute(attr, "BORDERINVALID", buffer);
      int borderinvalid = atoi(buffer);
      getAttribute(attr, "LARGEINVALID", buffer);
      int largeinvalid  = atoi(buffer);
      
      xmlparsedata->det->grading = 
	newGenEventGrading(invalid, borderinvalid, largeinvalid,
			   &xmlparsedata->status);

    } else if (!strcmp(Uelement, "GRADE")) {
      char buffer[MAXMSG]; // String buffer.

      getAttribute(attr, "P11", buffer);
      int p11  = atoi(buffer);
      getAttribute(attr, "P12", buffer);
      int p12  = atoi(buffer);
      getAttribute(attr, "P13", buffer);
      int p13  = atoi(buffer);

      getAttribute(attr, "P21", buffer);
      int p21  = atoi(buffer);
      getAttribute(attr, "P23", buffer);
      int p23  = atoi(buffer);

      getAttribute(attr, "P31", buffer);
      int p31  = atoi(buffer);
      getAttribute(attr, "P32", buffer);
      int p32  = atoi(buffer);
      getAttribute(attr, "P33", buffer);
      int p33  = atoi(buffer);

      getAttribute(attr, "GRADE", buffer);
      int grade = atoi(buffer);
      
      GenEventGrade* ggrade = newGenEventGrade(p11, p12, p13,
					       p21, p23,
					       p31, p32, p33,
					       grade, 
					       &xmlparsedata->status);
      if (EXIT_SUCCESS!=xmlparsedata->status) return;
      addGenEventGrade(xmlparsedata->det->grading,
		       ggrade,
		       &xmlparsedata->status);
      if (EXIT_SUCCESS!=xmlparsedata->status) return;
    }

    else { // Elements with independent attributes.

      // Loop over the different attributes.
      int i;
      for (i=0; attr[i]; i+=2) {
      
	// Convert the attribute to an upper case string.
	strcpy(Uattribute, attr[i]);
	strtoupper(Uattribute);

	// Check the XML element name.
	if (!strcmp(Uelement, "DIMENSIONS")) {
	  if (!strcmp(Uattribute, "XWIDTH")) {
	    xmlparsedata->det->pixgrid->xwidth = atoi(attr[i+1]);
	  } else if (!strcmp(Uattribute, "YWIDTH")) {
	    xmlparsedata->det->pixgrid->ywidth = atoi(attr[i+1]);
	  }
	}
      
	else if (!strcmp(Uelement, "WCS")) {
	  if (!strcmp(Uattribute, "XRPIX")) {
	    xmlparsedata->det->pixgrid->xrpix = (float)atof(attr[i+1]);
	  } else if (!strcmp(Uattribute, "YRPIX")) {
	    xmlparsedata->det->pixgrid->yrpix = (float)atof(attr[i+1]);
	  } else if (!strcmp(Uattribute, "XRVAL")) {
	    xmlparsedata->det->pixgrid->xrval = (float)atof(attr[i+1]);
	  } else if (!strcmp(Uattribute, "YRVAL")) {
	    xmlparsedata->det->pixgrid->yrval = (float)atof(attr[i+1]);
	  } else if (!strcmp(Uattribute, "XDELT")) {
	    xmlparsedata->det->pixgrid->xdelt = (float)atof(attr[i+1]);
	  } else if (!strcmp(Uattribute, "YDELT")) {
	    xmlparsedata->det->pixgrid->ydelt = (float)atof(attr[i+1]);
	  }
	}
	
	else if (!strcmp(Uelement, "PIXELBORDER")) {
	  if (!strcmp(Uattribute, "X")) {
	    xmlparsedata->det->pixgrid->xborder = (float)atof(attr[i+1]);
	  } else if (!strcmp(Uattribute, "Y")) {
	    xmlparsedata->det->pixgrid->yborder = (float)atof(attr[i+1]);
	  }
	}

	else if (!strcmp(Uelement, "RMF")) {
	  if (!strcmp(Uattribute, "FILENAME")) {
	    // Load the detector response file (RSP/RMF).
	    char buffer[MAXMSG];
	    strcpy(buffer, attr[i+1]);
	    xmlparsedata->det->rmf = loadRMF(buffer, &xmlparsedata->status);
	  }
	}

	else if (!strcmp(Uelement, "ARF")) {
	  if (!strcmp(Uattribute, "FILENAME")) {
	    // Load the detector ARF.
	    char buffer[MAXMSG];
	    strcpy(buffer, attr[i+1]);
	    xmlparsedata->det->arf = loadARF(buffer, &xmlparsedata->status);
	  }
	}
      
	else if (!strcmp(Uelement, "PSF")) {
	  // The focal length must be specified before load the PSF.
	  // Check if this is the case.
	  if (xmlparsedata->det->focal_length<=0.) {
	    xmlparsedata->status=EXIT_FAILURE;
	    HD_ERROR_THROW("Error: Telescope focal length must be specified "
			   "before loading the PSF!\n", xmlparsedata->status);
	    return;
	  }
	  if (!strcmp(Uattribute, "FILENAME")) {
	    // Load the PSF.
	    char buffer[MAXMSG];
	    strcpy(buffer, attr[i+1]);
	    xmlparsedata->det->psf = newPSF(buffer,
					    xmlparsedata->det->focal_length,
					    &xmlparsedata->status);
	  }
	}

	else if (!strcmp(Uelement, "CODEDMASK")) {
	  if (!strcmp(Uattribute, "FILENAME")) {
	    // Load the CodedMask.
	    char buffer[MAXMSG];
	    strcpy(buffer, attr[i+1]);
	    xmlparsedata->det->coded_mask = getCodedMaskFromFile(buffer, &xmlparsedata->status);
	  }
	}

	else if (!strcmp(Uelement, "VIGNETTING")) {
	  if (!strcmp(Uattribute, "FILENAME")) {
	    // Load the Vignetting function.
	    char buffer[MAXMSG];
	    strcpy(buffer, attr[i+1]);
	    xmlparsedata->det->vignetting = newVignetting(buffer, &xmlparsedata->status);
	  }
	}

	else if (!strcmp(Uelement, "FOCALLENGTH")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmlparsedata->det->focal_length = (float)atof(attr[i+1]);
	  }
	}

	else if (!strcmp(Uelement, "FOV")) {
	  if (!strcmp(Uattribute, "DIAMETER")) {
	    xmlparsedata->det->fov_diameter = (float)(atof(attr[i+1])*M_PI/180.);
	  }
	}

	else if (!strcmp(Uelement, "CTE")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmlparsedata->det->cte = (float)atof(attr[i+1]);
	  }
	}
	
	else if (!strcmp(Uelement, "BADPIXMAP")) {
	  if (!strcmp(Uattribute, "FILENAME")) {
	    // Load the detector bad pixel map.
	    char buffer[MAXMSG];
	    strcpy(buffer, attr[i+1]);
	    xmlparsedata->det->badpixmap = loadBadPixMap(buffer, &xmlparsedata->status);
	  }
	}

	else if (!strcmp(Uelement, "SPLIT")) {
	  if (!strcmp(Uattribute, "TYPE")) {
	    strcpy(Uvalue, attr[i+1]);
	    strtoupper(Uvalue);
	    if (!strcmp(Uvalue, "NONE")) {
	      xmlparsedata->det->split->type = GS_NONE;
	    } else if (!strcmp(Uvalue, "GAUSS")) {
	      xmlparsedata->det->split->type = GS_GAUSS;
	    } else if (!strcmp(Uvalue, "EXPONENTIAL")) {
	      xmlparsedata->det->split->type = GS_EXPONENTIAL;
	    }
	  } else if (!strcmp(Uattribute, "PAR1")) {
	    xmlparsedata->det->split->par1 = atof(attr[i+1]);
	  }
	}

	else if (!strcmp(Uelement, "EVENTFILE")) {
	  if (!strcmp(Uattribute, "TEMPLATE")) {
	    strcpy(xmlparsedata->det->eventfile_template, attr[i+1]);
	  }
	}

	else if (!strcmp(Uelement, "PATTERNFILE")) {
	  if (!strcmp(Uattribute, "TEMPLATE")) {
	    strcpy(xmlparsedata->det->patternfile_template, attr[i+1]);
	  }
	}

	else if (!strcmp(Uelement, "READOUT")) {
	  if (!strcmp(Uattribute, "MODE")) {
	    strcpy(Uvalue, attr[i+1]);
	    strtoupper(Uvalue);
	    if (!strcmp(Uvalue, "TIME")) {
	      xmlparsedata->det->readout_trigger = GENDET_TIME_TRIGGERED;
	    } else if (!strcmp(Uvalue, "EVENT")) {
	      xmlparsedata->det->readout_trigger = GENDET_EVENT_TRIGGERED;
	    }
	  }
	}
      
	else if (!strcmp(Uelement, "WAIT")) {
	  if (!strcmp(Uattribute, "TIME")) {
	    CLWait* clwait = newCLWait(atof(attr[i+1]), &xmlparsedata->status);
	    append2ClockList(xmlparsedata->det->clocklist, CL_WAIT, 
			     clwait, &xmlparsedata->status);
	  }
	}
	
	else if (!strcmp(Uelement, "CLEARLINE")) {
	  if (!strcmp(Uattribute, "LINEINDEX")) {
	    CLClearLine* clclearline = newCLClearLine(atoi(attr[i+1]), 
						      &xmlparsedata->status);
	    append2ClockList(xmlparsedata->det->clocklist, CL_CLEARLINE, 
			     clclearline, &xmlparsedata->status);
	  }
	}

	else if (!strcmp(Uelement, "THRESHOLD_READOUT_LO_KEV")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmlparsedata->det->threshold_readout_lo_keV = (float)atof(attr[i+1]);
	    headas_chat(3, "lower readout threshold: %.3lf keV\n", 
			xmlparsedata->det->threshold_readout_lo_keV);
	  }
	}
	
	else if (!strcmp(Uelement, "THRESHOLD_READOUT_UP_KEV")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmlparsedata->det->threshold_readout_up_keV = (float)atof(attr[i+1]);
	    headas_chat(3, "upper readout threshold: %.3lf keV\n", 
			xmlparsedata->det->threshold_readout_up_keV);
	  }
	}
	
	else if (!strcmp(Uelement, "THRESHOLD_READOUT_LO_PHA")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmlparsedata->det->threshold_readout_lo_PHA = (long)atoi(attr[i+1]);
	    headas_chat(3, "lower readout threshold: %ld PHA\n", 
			xmlparsedata->det->threshold_readout_lo_PHA);
	  }
	}
	
	else if (!strcmp(Uelement, "THRESHOLD_READOUT_UP_PHA")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmlparsedata->det->threshold_readout_up_PHA = (long)atoi(attr[i+1]);
	    headas_chat(3, "upper readout threshold: %ld PHA\n", 
			xmlparsedata->det->threshold_readout_up_PHA);
	  }
	}

	else if (!strcmp(Uelement, "THRESHOLD_EVENT_LO_KEV")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmlparsedata->det->threshold_event_lo_keV = (float)atof(attr[i+1]);
	    headas_chat(3, "lower event threshold: %.3lf keV\n", 
			xmlparsedata->det->threshold_event_lo_keV);
	  }
	}

	else if (!strcmp(Uelement, "THRESHOLD_SPLIT_LO_KEV")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmlparsedata->det->threshold_split_lo_keV = (float)atof(attr[i+1]);
	    headas_chat(3, "lower split threshold: %.3lf keV\n", 
			xmlparsedata->det->threshold_split_lo_keV);
	  }
	}

	else if (!strcmp(Uelement, "THRESHOLD_SPLIT_LO_FRACTION")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmlparsedata->det->threshold_split_lo_fraction = (float)atof(attr[i+1]);
	    headas_chat(3, "lower split threshold: %.1lf %%\n", 
			xmlparsedata->det->threshold_split_lo_fraction*100.);
	  }
	}

	if (EXIT_SUCCESS!=xmlparsedata->status) return;
      } 
      // END of loop over different attributes.
    }
    // END of elements with independent attributes
  }
  // END of elements with attributes.
}



static void GenDetXMLElementEnd(void* parsedata, const char* el) 
{
  struct XMLParseData* xmlparsedata = (struct XMLParseData*)parsedata;

  (void)el; // Unused parameter.

  // Check if an error has occurred previously.
  if (EXIT_SUCCESS!=xmlparsedata->status) return;

  return;
}



int addGenDetPhotonImpact(GenDet* const det, 
			  const Impact* const impact, 
			  int* const status)
{
  // Call the detector operating clock routine.
  operateGenDetClock(det, impact->time, status);
  if (EXIT_SUCCESS!=*status) return(0);

  // Determine the measured detector channel (PHA channel) according 
  // to the RMF.
  // The channel is obtained from the RMF using the corresponding
  // HEAdas routine which is based on drawing a random number.
  long channel;
  ReturnChannel(det->rmf, impact->energy, 1, &channel);

  // Check if the photon is really measured. If the
  // PHA channel returned by the HEAdas RMF function is '-1', 
  // the photon is not detected.
  // This can happen, if the RMF actually is an RSP, i.e. it 
  // includes ARF contributions, e.g., 
  // the detector quantum efficiency and filter transmission.
  if (0>channel) {
    return(0); // Break the function (photon is not detected).
  }

  // Get the corresponding created charge.
  // NOTE: In this simulation the charge is represented by the nominal
  // photon energy [keV] which corresponds to the PHA channel according 
  // to the EBOUNDS table.
  float charge = getEBOUNDSEnergy(channel, det->rmf, 0);
  assert(charge>=0.);

  // Create split events.
  int npixels = makeGenSplitEvents(det->split, &impact->position, charge, 
				   det->pixgrid, det->line)>0;
  
  // Return the number of affected pixels.
  return(npixels);
}



void operateGenDetClock(GenDet* const det, const double time, int* const status)
{
  // Check if the detector operation setup is time-triggered.
  if (GENDET_TIME_TRIGGERED!=det->readout_trigger) return;

  // Get the next element from the clock list.
  CLType type;
  void* element=NULL;
  do {
    CLReadoutLine* clreadoutline=NULL;
    CLClearLine*   clclearline  =NULL;
    CLWait* clwait              =NULL;

    getClockListElement(det->clocklist, time, &type, &element, status);
    if (EXIT_SUCCESS!=*status) return;

    switch (type) {
    case CL_NONE:
      // No operation has to be performed. The clock list is
      // currently in a wait status.
      break;
    case CL_WAIT:
      // A waiting period is finished.
      clwait = (CLWait*)element;
      // Apply the bad pixel map (if available) with the bad pixel values weighted 
      // by the waiting time.
      if (NULL!=det->badpixmap) {
	applyBadPixMap(det->badpixmap, clwait->time, encounterGenDetBadPix, det->line);
      }
      break;
    case CL_LINESHIFT:
      GenDetLineShift(det);
      break;
    case CL_READOUTLINE:
      clreadoutline = (CLReadoutLine*)element;
      GenDetReadoutLine(det, clreadoutline->lineindex, 
			clreadoutline->readoutindex, 
			status);
      break;
    case CL_CLEARLINE:
      clclearline = (CLClearLine*)element;
      GenDetClearLine(det, clclearline->lineindex);
      break;
    }
    if(EXIT_SUCCESS!=*status) return;
  } while(type!=CL_NONE);
}



void GenDetLineShift(GenDet* const det)
{
  headas_chat(5, "lineshift\n");

  // Check if the detector contains more than 1 line.
  if (2>det->pixgrid->ywidth) return;

  // Apply the Charge Transfer Efficiency.
  int ii;
  if (det->cte!=1.) {
#pragma omp for 
    for (ii=1; ii<det->pixgrid->ywidth; ii++) {
      if (0!=det->line[ii]->anycharge) {
	int jj;
	for (jj=0; jj<det->line[ii]->xwidth; jj++) {
	  det->line[ii]->charge[jj] *= det->cte;
	  // TODO RM
	  printf("ii=%d jj=%d\n", ii, jj);
	}
      }
    }
  }

  // Add the charges in line 1 to line 0.
  addGenDetLine(det->line[0], det->line[1]);

  // Clear the charges in line 1, as they are now contained in line 0.
  clearGenDetLine(det->line[1]);

  // Shift the other lines in increasing order and put the newly cleared 
  // original line number 1 at the end as the last line.
  GenDetLine* buffer = det->line[1];
  for (ii=1; ii<det->pixgrid->ywidth-1; ii++) {
    det->line[ii] = det->line[ii+1];
  }
  det->line[det->pixgrid->ywidth-1] = buffer;
}



void GenDetReadoutLine(GenDet* const det, const int lineindex, 
		       const int readoutindex, int* const status)
{
  headas_chat(5, "read out line %d as %d\n", lineindex, readoutindex);

  // Event data structure.
  GenEvent event = {.time=0.};

  while (readoutGenDetLine(det->line[lineindex], &event)) {

    // Apply the charge thresholds.
    if (event.charge<=det->threshold_readout_lo_keV) continue;

    if (det->threshold_readout_up_keV >= 0.) {
      if (event.charge>=det->threshold_readout_up_keV) continue;
    }

    // Apply the detector response.
    event.pha = getEBOUNDSChannel(event.charge, det->rmf);

    // Store the additional information.
    event.rawy  = readoutindex;
    event.time  = det->clocklist->time;  // Time of detection.
    event.frame = det->clocklist->frame; // Frame of detection.

    // Store the event in the output event file.
    addGenEvent2File(det->eventfile, &event, status);
  }
}



void GenDetClearLine(GenDet* const det, const int lineindex) {
  clearGenDetLine(det->line[lineindex]);
}



void GenDetNewEventFile(GenDet* const det, const char* const filename, 
			int* const status)
{
  // Check if there already is an open event file.
  if (NULL!=det->eventfile) {
    // Close the file
    destroyGenEventFile(&det->eventfile, status);
  }

  // Filename of the template file.
  char template[MAXMSG];
  // Get the name of the FITS template directory.
  // First try to read it from the environment variable.
  // If the variable does not exist, read it from the PIL.
  char* buffer;
  if (NULL!=(buffer=getenv("SIXT_FITS_TEMPLATES"))) {
    strcpy(template, buffer);
  } else {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not read environment variable 'SIXT_FITS_TEMPLATES'!\n", 
		   *status);
    return;
  }

  // Append the filename of the template file itself.
  strcat(template, "/");
  strcat(template, det->eventfile_template);

  // Open a new event file from the specified template.
  det->eventfile = openNewGenEventFile(filename, template, status);

  // Insert the header keywords describing the detector.
  if (fits_update_key(det->eventfile->fptr, TINT, "NXDIM", 
		      &det->pixgrid->xwidth, "number of pixels in x-direction", 
		      status)) return;
  if (fits_update_key(det->eventfile->fptr, TINT, "NYDIM", 
		      &det->pixgrid->ywidth, "number of pixels in y-direction", 
		      status)) return;
  if (fits_update_key(det->eventfile->fptr, TLONG, "DETCHANS",
		      &det->rmf->NumberChannels, "number of EBOUNDS channels",
		      status)) return;
  if (fits_update_key(det->eventfile->fptr, TLONG, "TLMIN1",
  		      &det->rmf->FirstChannel, "first channel of EBOUNDS",
  		      status)) return;
  long tlmax1 = det->rmf->FirstChannel + det->rmf->NumberChannels - 1;
  if (fits_update_key(det->eventfile->fptr, TLONG, "TLMAX1",
  		      &tlmax1, "last channel of EBOUNDS",
  		      status)) return;
}



static void expandXML(struct XMLBuffer* const buffer, int* const status)
{
  struct XMLPreParseData data;
  
  do {

    // Parse XML code in the xmlbuffer using the expat library.
    // Get an XML_Parser object.
    XML_Parser parser = XML_ParserCreate(NULL);
    if (NULL==parser) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: Could not allocate memory for XML parser!\n", *status);
      return;
    }

    // Set data that is passed to the handler functions.
    XML_SetUserData(parser, &data);

    // Set the handler functions.
    XML_SetElementHandler(parser, expandXMLElementStart, expandXMLElementEnd);

    // Set initial values.
    data.further_loops = 0;
    data.loop_depth = 0;
    data.loop_start = 0;
    data.loop_end   = 0;
    data.loop_increment = 0;
    data.output_buffer = newXMLBuffer(status);
    data.loop_buffer   = newXMLBuffer(status);
    data.status = EXIT_SUCCESS;

    // Process all the data in the string buffer.
    const int done=1;
    if (!XML_Parse(parser, buffer->text, strlen(buffer->text), done)) {
      // Parse error.
      *status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "Error: Parsing XML code failed:\n%s\n", 
	      XML_ErrorString(XML_GetErrorCode(parser)));
      printf("%s", buffer->text);
      HD_ERROR_THROW(msg, *status);
      return;
    }
    // Check for errors.
    if (EXIT_SUCCESS!=data.status) {
      *status = data.status;
      return;
    }

    // Copy the output XMLBuffer to the input XMLBuffer ...
    copyXMLBuffer(buffer, data.output_buffer, status);
    if (EXIT_SUCCESS!=*status) return;
    // ... and release allocated memory.
    destroyXMLBuffer(&data.output_buffer);
    destroyXMLBuffer(&data.loop_buffer);

    XML_ParserFree(parser);

  } while (data.further_loops>0);

  // Replace arithmetic +/- expressions.
  execArithmeticOpsInXMLBuffer(buffer, status);
}



static void expandXMLElementStart(void* data, const char* el, 
				  const char** attr)
{
  struct XMLPreParseData* mydata = (struct XMLPreParseData*)data;

  // Pointer to the right output buffer (either mydata->output_buffer 
  // or mydata->loop_buffer).
  struct XMLBuffer* output=mydata->output_buffer;

  // Convert the element to an upper case string.
  char Uelement[MAXMSG]; // Upper case version of XML element
  strcpy(Uelement, el);
  strtoupper(Uelement);

  // Check if the element is a loop tag.
  if (!strcmp(Uelement, "LOOP")) {

    // Check if this is the outermost loop.
    if (0==mydata->loop_depth) {
      // Read the loop parameters.
      mydata->loop_start    =0;
      mydata->loop_end      =0;
      mydata->loop_increment=0;
      mydata->loop_variable[0]='\0';

      int ii=0;
      while (attr[ii]) {
	char Uattribute[MAXMSG];
	char Uvalue[MAXMSG];
	strcpy(Uattribute, attr[ii]);
	strtoupper(Uattribute);
	strcpy(Uvalue, attr[ii+1]);
	strtoupper(Uvalue);

	if (!strcmp(Uattribute, "START")) {
	  mydata->loop_start = atoi(attr[ii+1]);
	} else if (!strcmp(Uattribute, "END")) {
	  mydata->loop_end = atoi(attr[ii+1]);
	} else if (!strcmp(Uattribute, "INCREMENT")) {
	  mydata->loop_increment = atoi(attr[ii+1]);
	} else if (!strcmp(Uattribute, "VARIABLE")) {
	  strcpy(mydata->loop_variable, attr[ii+1]);
	}

	ii+=2;
      }
      // END of loop over all attributes.
      
      // Check if parameters are set to valid values.
      if ((mydata->loop_end-mydata->loop_start)*mydata->loop_increment<=0) {
	mydata->status = EXIT_FAILURE;
	HD_ERROR_THROW("Error: Invalid XML loop parameters!\n", mydata->status);
	return;
      }
      mydata->loop_depth++;
      return;

    } else {
      // Inner loop.
      mydata->further_loops = 1;
      mydata->loop_depth++;
    }
    // END of check if this is the outermost loop or an inner loop.
  }
  // END of check for loop tag.

  
  // If we are inside a loop, print to the loop buffer.
  if (mydata->loop_depth>0) {
    output=mydata->loop_buffer;
  }

  // Print the start tag to the right buffer.
  char buffer[MAXMSG];
  if (sprintf(buffer, "<%s", el) >= MAXMSG) {
    mydata->status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: XML element string too long!\n", EXIT_FAILURE);
    return;
  }
  addString2XMLBuffer(output, buffer, &mydata->status);
  if (EXIT_SUCCESS!=mydata->status) return;

  int ii=0;
  while(attr[ii]) {
    if (sprintf(buffer, " %s=\"%s\"", attr[ii], attr[ii+1]) >= MAXMSG) {
      mydata->status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: XML element string too long!\n", EXIT_FAILURE);
      return;
    }
    addString2XMLBuffer(output, buffer, &mydata->status);
    if (EXIT_SUCCESS!=mydata->status) return;
    
    ii+=2;
  }

  addString2XMLBuffer(output, ">", &mydata->status);
  if (EXIT_SUCCESS!=mydata->status) return;
}



static void expandXMLElementEnd(void* data, const char* el)
{
  struct XMLPreParseData* mydata = (struct XMLPreParseData*)data;

  // Pointer to the right output buffer (either mydata->output_buffer 
  // or mydata->loop_buffer).
  struct XMLBuffer* output=mydata->output_buffer;

  // Convert the element to an upper case string.
  char Uelement[MAXMSG]; // Upper case version of XML element
  strcpy(Uelement, el);
  strtoupper(Uelement);

  // Check if the element is a loop end tag.
  if (!strcmp(Uelement, "LOOP")) {

    // Check if the outer loop is finished.
    // In that case add the loop buffer n-times to the output buffer.
    if (1==mydata->loop_depth) {
      int ii;
      for (ii=mydata->loop_start; 
	   ((ii<=mydata->loop_end)&&(mydata->loop_increment>0)) ||
	     ((ii>=mydata->loop_end)&&(mydata->loop_increment<0)); 
	   ii+=mydata->loop_increment) {
	// Copy loop buffer to separate XMLBuffer before replacing the
	// variables, since they have to be preserved for the following loop
	// repetitions.
	struct XMLBuffer* replacedBuffer=newXMLBuffer(&mydata->status);
	if (EXIT_SUCCESS!=mydata->status) return;
	copyXMLBuffer(replacedBuffer, mydata->loop_buffer, &mydata->status);
	if (EXIT_SUCCESS!=mydata->status) return;

	// Replace $variables by integer values.
	if (strlen(replacedBuffer->text)>0) {
	  char stringvalue[MAXMSG];
	  sprintf(stringvalue, "%d", ii);
	  replaceInXMLBuffer(replacedBuffer, mydata->loop_variable,
			     stringvalue, &mydata->status);
	  if (EXIT_SUCCESS!=mydata->status) return;
	}

	// Add the loop content to the output buffer.
	addString2XMLBuffer(mydata->output_buffer, replacedBuffer->text, 
			    &mydata->status);
	if (EXIT_SUCCESS!=mydata->status) return;

	destroyXMLBuffer(&replacedBuffer);
      }
      // Clear the loop buffer.
      destroyXMLBuffer(&mydata->loop_buffer);
      mydata->loop_buffer=newXMLBuffer(&mydata->status);
      if (EXIT_SUCCESS!=mydata->status) return;
      
      // Now we are outside of any loop.
      mydata->loop_depth--;
      return;

    } else {
      // We are still inside some outer loop.
      mydata->loop_depth--;
    }
  }
  
  // If we are inside a loop, print to the loop buffer.
  if (mydata->loop_depth>0) {
    output=mydata->loop_buffer;
  }

  // Print the end tag to the right buffer.
  char buffer[MAXMSG];
  if (sprintf(buffer, "</%s>", el) >= MAXMSG) {
    mydata->status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: XML string element too long!\n", EXIT_FAILURE);
    return;
  }
  addString2XMLBuffer(output, buffer, &mydata->status);
  if (EXIT_SUCCESS!=mydata->status) return;
}



static void replaceInXMLBuffer(struct XMLBuffer* const buffer, 
			       const char* const old,
			       const char* const new, 
			       int* const status)
{
  // Pointer to the first occurence of the old string in the buffer text.
  char* occurence; 
  // Length of the old string.
  int len_old = strlen(old);

  if ((0==strlen(old)) || (0==strlen(buffer->text))) return;

  while (NULL!=(occurence=strstr(buffer->text, old))) {
    // Get the length of the tail.
    int len_tail = strlen(occurence)-len_old;

    // String tail after the first occurence of the old string in the buffer text.
    char* tail=(char*)malloc((1+len_tail)*sizeof(char));
    if (NULL==tail) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: Memory allocation for string buffer in XML "
		     "pre-parsing failed!\n", *status);
      return;
    }
    // Copy the tail of the string without the old string to the buffer.
    strcpy(tail, &occurence[len_old]);
    // Truncate the buffer string directly before the occurence of the old string.
    occurence[0] = '\0';
    // Append the new string to the truncated buffer string.
    addString2XMLBuffer(buffer, new, status);
    // Append the tail after the newly inserted new string.
    addString2XMLBuffer(buffer, tail, status);

    // Release memory of the tail buffer.
    free(tail);
  }
}



static void execArithmeticOpsInXMLBuffer(struct XMLBuffer* const buffer,
					 int* const status)
{
  char* occurrence=buffer->text;

  // Loop while a "+" or "-" sign is found in the buffer text.
  while (NULL!=(occurrence=strpbrk(occurrence, "+-"))) {

    // 1. Determine the first term.
    char* start = occurrence-1;
    // Scan forward until reaching a non-digit.
    while (strpbrk(start, "0123456789")==start) {
      start--;
    }
    start++;

    // Check if there really is a numeric term in front of the "+" or "-" sign.
    // If not (e.g. "e-4") continue with the next loop iteration.
    if (start==occurrence) {
      occurrence++;
      continue;
    }

    // Store the first value in a separate string.
    char svalue[MAXMSG];
    int ii;
    for (ii=0; start+ii<occurrence; ii++) {
      svalue[ii] = start[ii];
    }
    svalue[ii] = '\0';
    
    // Convert the string to an integer value.
    int ivalue1 = atoi(svalue);


    // 2. Determine the second term.
    char* end = occurrence+1;
    // Scan backward until reaching a non-digit.
    while (strpbrk(end, "0123456789")==end) {
      end++;
    }
    end--;

    // Check if there really is a numeric term behind of the "+" or "-" sign.
    // If not continue with the next loop iteration.
    if (end==occurrence) {
      occurrence++;
      continue;
    }

    // Store the second value in a separate string.
    for (ii=0; occurrence+ii<end; ii++) {
      svalue[ii] = occurrence[1+ii];
    }
    svalue[ii] = '\0';

    // Convert the string to an integer value.
    int ivalue2 = atoi(svalue);


    // Perform the arithmetic operation.
    int result=0;
    if (occurrence[0]=='+') {
      result = ivalue1 + ivalue2;
    } else if (occurrence[0]=='-') {
      result = ivalue1 - ivalue2;
    }


    // Store the result at the right position in the XMLBuffer text.
    // Determine the new sub-string.
    sprintf(svalue, "%d", result);

    // Get the length of the tail.
    int len_tail = strlen(end)-1;

    // String tail after the first occurrence of the old string in the buffer text.
    char* tail=(char*)malloc((1+len_tail)*sizeof(char));
    if (NULL==tail) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: Memory allocation for string buffer in XML "
		     "pre-parsing failed!\n", *status);
      return;
    }
    // Copy the tail of the string without the old string to the buffer.
    strcpy(tail, &end[1]);
    // Truncate the buffer string directly before the occurrence of the old string.
    start[0] = '\0';
    // Append the new string to the truncated buffer string.
    addString2XMLBuffer(buffer, svalue, status);
    // Append the tail after the newly inserted new string.
    addString2XMLBuffer(buffer, tail, status);

    // Release memory of the tail buffer.
    free(tail);
  }
  // END of loop over all occurrences of "+" or "-" signs.
}



void encounterGenDetBadPix(void* const data, const int x, const int y, const float value) {
  // Array of pointers to pixel lines.
  GenDetLine** line = (GenDetLine**)data;

  // Check if the bad pixel type.
  if (value < 0.) { // The pixel is a cold one.
    // Delete the charge in the pixel.
    line[y]->charge[x] = 0.;
  } else if (value > 0.) { // The pixel is a hot one.
    // Add additional charge to the pixel.
    addGenDetCharge2Pixel(line[y], x, value);
  }
}

