#include "gendet.h"


////////////////////////////////////////////////////////////////////
// Local data type declarations.
////////////////////////////////////////////////////////////////////


/** Data structure given to the XML handler to transfer data. */
struct XMLParseData {
  GenDet* det;
  int status;
};


////////////////////////////////////////////////////////////////////
// Static function declarations
////////////////////////////////////////////////////////////////////


/** Handler for the start of an XML element. */
static void GenDetXMLElementStart(void* data, const char* el, 
				  const char** attr);
/** Handler for the end of an XML element. */
static void GenDetXMLElementEnd(void* data, const char* el);


////////////////////////////////////////////////////////////////////
// Program Code.
////////////////////////////////////////////////////////////////////


void parseGenDetXML(GenDet* const det, 
		    const char* const filename, 
		    int* const status)
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


  // Release memory.
  XML_ParserFree(parser);

  // Remove the XML string buffer.
  freeXMLBuffer(&xmlbuffer);


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
    headas_printf("*** warning: no specification for response file (RMF/RSP) "
		  "in XML definition ***");
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
  if (NULL!=det->rmf) {
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
  CHECK_STATUS_VOID(xmlparsedata->status);

  // Convert the element to an upper case string.
  strcpy(Uelement, el);
  strtoupper(Uelement);

  // Elements without attributes.
  if (!strcmp(Uelement, "LINESHIFT")) {
    CLLineShift* cllineshift = newCLLineShift(&xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);
    append2ClockList(xmlparsedata->det->clocklist, CL_LINESHIFT, 
		     cllineshift, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

  } else if (!strcmp(Uelement, "NEWFRAME")) {
    CLNewFrame* clnewframe = newCLNewFrame(&xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);
    append2ClockList(xmlparsedata->det->clocklist, CL_NEWFRAME, 
		     clnewframe, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

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
	    char buffer[MAXFILENAME];
	    strcpy(buffer, xmlparsedata->det->filepath);
	    strcat(buffer, attr[i+1]);
	    xmlparsedata->det->rmf = loadRMF(buffer, &xmlparsedata->status);
	  }
	}

	else if (!strcmp(Uelement, "ARF")) {
	  if (!strcmp(Uattribute, "FILENAME")) {
	    // Load the detector ARF.
	    char buffer[MAXFILENAME];
	    strcpy(buffer, xmlparsedata->det->filepath);
	    strcat(buffer, attr[i+1]);
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
	    char buffer[MAXFILENAME];
	    strcpy(buffer, xmlparsedata->det->filepath);
	    strcat(buffer, attr[i+1]);
	    xmlparsedata->det->psf = newPSF(buffer,
					    xmlparsedata->det->focal_length,
					    &xmlparsedata->status);
	  }
	}

	else if (!strcmp(Uelement, "CODEDMASK")) {
	  if (!strcmp(Uattribute, "FILENAME")) {
	    // Load the CodedMask.
	    char buffer[MAXFILENAME];
	    strcpy(buffer, xmlparsedata->det->filepath);
	    strcat(buffer, attr[i+1]);
	    xmlparsedata->det->coded_mask = 
	      getCodedMaskFromFile(buffer, &xmlparsedata->status);
	  }
	}

	else if (!strcmp(Uelement, "VIGNETTING")) {
	  if (!strcmp(Uattribute, "FILENAME")) {
	    // Load the Vignetting function.
	    char buffer[MAXFILENAME];
	    strcpy(buffer, xmlparsedata->det->filepath);
	    strcat(buffer, attr[i+1]);
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
	    char buffer[MAXFILENAME];
	    strcpy(buffer, xmlparsedata->det->filepath);
	    strcat(buffer, attr[i+1]);
	    xmlparsedata->det->badpixmap = loadBadPixMap(buffer, &xmlparsedata->status);
	  }
	}

	else if (!strcmp(Uelement, "EROBACKGROUND")) {
	  if (!strcmp(Uattribute, "FILENAME")) {
	    // Load the detector background model for
	    // cosmic rays.
	    char buffer[MAXFILENAME];
	    strcpy(buffer, xmlparsedata->det->filepath);
	    strcat(buffer, attr[i+1]);
	    eroBkgInitialize(buffer, &xmlparsedata->status);
	    xmlparsedata->det->erobackground = 1;
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


