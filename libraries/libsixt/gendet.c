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
  int status;
};


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
  det->rmf =NULL;
  det->psf =NULL;
  det->vignetting=NULL;
  det->clocklist =NULL;
  det->eventfile =NULL;

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

    // Close the event file.
    destroyGenEventFile(&(*det)->eventfile, status);

    // Free the PSF.
    free_psf(&(*det)->psf);
    
    // Free the vignetting Function.
    free_Vignetting(&(*det)->vignetting);

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

  // Set initial values before parsing the parameters from the XML file.
  det->pixgrid->xwidth =-1;
  det->pixgrid->ywidth =-1;
  det->pixgrid->xrpix  =-1.;
  det->pixgrid->yrpix  =-1.;
  det->pixgrid->xrval  =-1.;
  det->pixgrid->yrval  =-1.;
  det->pixgrid->xdelt  =-1.;
  det->pixgrid->ydelt  =-1.;
  det->readout_trigger = 0;
  det->cte             = 1.;
  det->threshold_readout_lo_PHA    = -1;
  det->threshold_readout_up_PHA    = -1;
  det->threshold_readout_lo_keV    =  0.;
  det->threshold_readout_up_keV    = -1.;
  det->threshold_event_lo_keV      = -1.;
  det->threshold_split_lo_keV      =  0.;
  det->threshold_split_lo_fraction =  0.;
  det->fov_diameter = 0.;
  det->focal_length = 0.;
  // Set string variables to empty strings.
  strcpy(det->eventfile_template, "");


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
    .det   = det,
    .status = EXIT_SUCCESS
  };
  XML_SetUserData(parser, &xmldata);

  // Set the handler functions.
  XML_SetElementHandler(parser, GenDetXMLElementStart, GenDetXMLElementEnd);

  int done=0;
  int len; // Number of chars in buffer.
  const int buffer_size=200;
  // Input buffer with an additional byte at the end for the 
  // termination of the string.
  char buffer[buffer_size+1];

  do {
    // Get a piece of input into the buffer.
    len = fread(buffer, 1, buffer_size, xmlfile);
    buffer[len]='\0'; // Terminate the string.
    done = feof(xmlfile);
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
    // Check for errors.
    if (EXIT_SUCCESS!=xmldata.status) {
      *status = xmldata.status;
      return;
    }
  } while (!done);
  XML_ParserFree(parser);

  // Check if all required parameters have been read successfully from 
  // the XML file.
  if (-1==det->pixgrid->xwidth) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for x-width of GenDet pixel array!\n", 
		   *status);
    return;    
  }  
  if (-1==det->pixgrid->ywidth) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for y-width of GenDet pixel array!\n", 
		   *status);
    return;    
  }

  if (0>det->pixgrid->xrpix) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for x reference pixel of GenDet!\n", 
		   *status);
    return;    
  }
  if (0>det->pixgrid->yrpix) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for y reference pixel of GenDet!\n", 
		   *status);
    return;    
  }

  if (0>det->pixgrid->xrval) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for x reference value of GenDet!\n", 
		   *status);
    return;    
  }
  if (0>det->pixgrid->yrval) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for y reference value of GenDet!\n", 
		   *status);
    return;    
  }

  if (0>det->pixgrid->xdelt) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for x pixel width of GenDet!\n", 
		   *status);
    return;    
  }
  if (0>det->pixgrid->ydelt) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for y pixel width of GenDet!\n", 
		   *status);
    return;    
  }
  
  if (NULL==det->rmf) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for response file of GenDet!\n", 
		   *status);
    return;    
  }

  if (NULL==det->psf) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for PSF!\n", 
		   *status);
    return;    
  }
  if (NULL==det->vignetting) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for Vignetting!\n", 
		   *status);
    return;    
  }
  if (0.>det->focal_length) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for the focal length of the telescope!\n", 
		   *status);
    return;    
  }
  if (0.>det->fov_diameter) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for the diameter of the telescope FoV!\n", 
		   *status);
    return;    
  }

  if (0==det->readout_trigger) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No specification found for the readout trigger of GenDet!\n", 
		   *status);
    return;
  }

  if (0==strlen(det->eventfile_template)) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No event file template specified!\n", *status);
    return;    
  }

  if (GS_NONE!=det->split->type) {
    if (det->split->par1<=0.) {
      *status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: No valid split model parameters!\n", *status);
      return;    
    }
  }

  if (0>det->threshold_event_lo_keV) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: No lower event threshold specified!\n", *status);
    return;
  }

  // END of checking, if all detector parameters have successfully been 
  // read from the XML file.

  // Close the file handler to the XML file.
  fclose(xmlfile);

  // If any thresholds have been specified in terms of PHA value,
  // set the corresponding charge threshold to the [keV] values
  // according to the RMF. If a charge threshold is given in addition,
  // its value is overwritten by the charge corresponding to the PHA 
  // specification. I.e., the PHA thresholds have a higher priority.
  if (det->threshold_readout_lo_PHA>-1) {
    det->threshold_readout_lo_keV = getEnergy(det->threshold_readout_lo_PHA, det->rmf, -1);
    headas_chat(3, "set lower readout threshold to %.3lf keV (PHA %ld)\n", 
		det->threshold_readout_lo_keV, det->threshold_readout_lo_PHA);
  }
  if (det->threshold_readout_up_PHA>-1) {
    det->threshold_readout_up_keV = getEnergy(det->threshold_readout_up_PHA, det->rmf,  1);
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



static void GenDetXMLElementStart(void* data, const char* el, const char** attr) 
{
  struct XMLData* xmldata = (struct XMLData*)data;
  char Uelement[MAXMSG];   // Upper case version of XML element
  char Uattribute[MAXMSG]; // Upper case version of XML attribute
  char Uvalue[MAXMSG];     // Upper case version of XML attribute value

  // Check if an error has occurred previously.
  if (EXIT_SUCCESS!=xmldata->status) return;

  // Convert the element to an upper case string.
  strcpy(Uelement, el);
  strtoupper(Uelement);

  // Elements without attributes.
  if (!strcmp(Uelement, "LINESHIFT")) {
    CLLineShift* cllineshift = newCLLineShift(&xmldata->status);
    append2ClockList(xmldata->det->clocklist, CL_LINESHIFT, 
		     cllineshift, &xmldata->status);
    
  } else { 
    
    // Elements with attributes.

    if (!strcmp(Uelement, "READOUTLINE")) {
      char buffer[MAXMSG]; // String buffer.
      getAttribute(attr, "LINEINDEX", buffer);
      int lineindex    = atoi(buffer);
      if (lineindex<0) {
	xmldata->status=EXIT_FAILURE;
	HD_ERROR_THROW("Error: Negative index for readout line!\n", xmldata->status);
	return;
      }
      getAttribute(attr, "READOUTINDEX", buffer);
      int readoutindex = atoi(buffer);
      if (readoutindex<0) {
	xmldata->status=EXIT_FAILURE;
	HD_ERROR_THROW("Error: Negative index for readout line!\n", xmldata->status);
	return;
      }
      CLReadoutLine* clreadoutline = newCLReadoutLine(lineindex,
						      readoutindex,
						      &xmldata->status);
      append2ClockList(xmldata->det->clocklist, CL_READOUTLINE, 
		       clreadoutline, &xmldata->status);
	
    } else { // Elements with independent attributes.

      // Loop over the different attributes.
      int i;
      for (i=0; attr[i]; i+=2) {
      
	// Convert the attribute to an upper case string.
	strcpy(Uattribute, attr[i]);
	strtoupper(Uattribute);

	// Check the XML element name.
	if (!strcmp(Uelement, "DIMENSIONS")) {
	  if (!strcmp(Uattribute, "XWIDTH")) {
	    xmldata->det->pixgrid->xwidth = atoi(attr[i+1]);
	  } else if (!strcmp(Uattribute, "YWIDTH")) {
	    xmldata->det->pixgrid->ywidth = atoi(attr[i+1]);
	  }
	}
      
	else if (!strcmp(Uelement, "WCS")) {
	  if (!strcmp(Uattribute, "XRPIX")) {
	    xmldata->det->pixgrid->xrpix = (float)atof(attr[i+1]);
	  } else if (!strcmp(Uattribute, "YRPIX")) {
	    xmldata->det->pixgrid->yrpix = (float)atof(attr[i+1]);
	  } else if (!strcmp(Uattribute, "XRVAL")) {
	    xmldata->det->pixgrid->xrval = (float)atof(attr[i+1]);
	  } else if (!strcmp(Uattribute, "YRVAL")) {
	    xmldata->det->pixgrid->yrval = (float)atof(attr[i+1]);
	  } else if (!strcmp(Uattribute, "XDELT")) {
	    xmldata->det->pixgrid->xdelt = (float)atof(attr[i+1]);
	  } else if (!strcmp(Uattribute, "YDELT")) {
	    xmldata->det->pixgrid->ydelt = (float)atof(attr[i+1]);
	  }
	}
	
	else if (!strcmp(Uelement, "RESPONSE")) {
	  if (!strcmp(Uattribute, "FILENAME")) {
	    // Load the detector response file.
	    char buffer[MAXMSG];
	    strcpy(buffer, attr[i+1]);
	    xmldata->det->rmf = loadRMF(buffer, &xmldata->status);
	  }
	}
      
	else if (!strcmp(Uelement, "PSF")) {
	  if (!strcmp(Uattribute, "FILENAME")) {
	    // Load the PSF.
	    char buffer[MAXMSG];
	    strcpy(buffer, attr[i+1]);
	    xmldata->det->psf = get_psf(buffer, &xmldata->status);
	  }
	}

	else if (!strcmp(Uelement, "VIGNETTING")) {
	  if (!strcmp(Uattribute, "FILENAME")) {
	    // Load the Vignetting function.
	    char buffer[MAXMSG];
	    strcpy(buffer, attr[i+1]);
	    xmldata->det->vignetting = get_Vignetting(buffer, &xmldata->status);
	  }
	}

	else if (!strcmp(Uelement, "FOCALLENGTH")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmldata->det->focal_length = (float)atof(attr[i+1]);
	  }
	}

	else if (!strcmp(Uelement, "FOV")) {
	  if (!strcmp(Uattribute, "DIAMETER")) {
	    xmldata->det->fov_diameter = (float)(atof(attr[i+1])*M_PI/180.);
	  }
	}

	else if (!strcmp(Uelement, "CTE")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmldata->det->cte = (float)atof(attr[i+1]);
	  }
	}
	
	else if (!strcmp(Uelement, "SPLIT")) {
	  if (!strcmp(Uattribute, "TYPE")) {
	    strcpy(Uvalue, attr[i+1]);
	    strtoupper(Uvalue);
	    if (!strcmp(Uvalue, "NONE")) {
	      xmldata->det->split->type = GS_NONE;
	    } else if (!strcmp(Uvalue, "GAUSS")) {
	      xmldata->det->split->type = GS_GAUSS;
	    } else if (!strcmp(Uvalue, "EXPONENTIAL")) {
	      xmldata->det->split->type = GS_EXPONENTIAL;
	    }
	  } else if (!strcmp(Uattribute, "PAR1")) {
	    xmldata->det->split->par1 = atof(attr[i+1]);
	  }
	}

	else if (!strcmp(Uelement, "EVENTFILE")) {
	  if (!strcmp(Uattribute, "TEMPLATE")) {
	    strcpy(xmldata->det->eventfile_template, attr[i+1]);
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
      
	else if (!strcmp(Uelement, "WAIT")) {
	  if (!strcmp(Uattribute, "TIME")) {
	    CLWait* clwait = newCLWait(atof(attr[i+1]), &xmldata->status);
	    append2ClockList(xmldata->det->clocklist, CL_WAIT, 
			     clwait, &xmldata->status);
	  }
	}
	
	else if (!strcmp(Uelement, "CLEARLINE")) {
	  if (!strcmp(Uattribute, "LINEINDEX")) {
	    CLClearLine* clclearline = newCLClearLine(atoi(attr[i+1]), 
						      &xmldata->status);
	    append2ClockList(xmldata->det->clocklist, CL_CLEARLINE, 
			     clclearline, &xmldata->status);
	  }
	}

	else if (!strcmp(Uelement, "THRESHOLD_READOUT_LO_KEV")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmldata->det->threshold_readout_lo_keV = (float)atof(attr[i+1]);
	    headas_chat(3, "lower readout threshold: %.3lf keV\n", 
			xmldata->det->threshold_readout_lo_keV);
	  }
	}
	
	else if (!strcmp(Uelement, "THRESHOLD_READOUT_UP_KEV")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmldata->det->threshold_readout_up_keV = (float)atof(attr[i+1]);
	    headas_chat(3, "upper readout threshold: %.3lf keV\n", 
			xmldata->det->threshold_readout_up_keV);
	  }
	}
	
	else if (!strcmp(Uelement, "THRESHOLD_READOUT_LO_PHA")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmldata->det->threshold_readout_lo_PHA = (long)atoi(attr[i+1]);
	    headas_chat(3, "lower readout threshold: %ld PHA\n", 
			xmldata->det->threshold_readout_lo_PHA);
	  }
	}
	
	else if (!strcmp(Uelement, "THRESHOLD_READOUT_UP_PHA")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmldata->det->threshold_readout_up_PHA = (long)atoi(attr[i+1]);
	    headas_chat(3, "upper readout threshold: %ld PHA\n", 
			xmldata->det->threshold_readout_up_PHA);
	  }
	}

	else if (!strcmp(Uelement, "THRESHOLD_EVENT_LO_KEV")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmldata->det->threshold_event_lo_keV = (float)atof(attr[i+1]);
	    headas_chat(3, "lower event threshold: %.3lf keV\n", 
			xmldata->det->threshold_event_lo_keV);
	  }
	}

	else if (!strcmp(Uelement, "THRESHOLD_SPLIT_LO_KEV")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmldata->det->threshold_split_lo_keV = (float)atof(attr[i+1]);
	    headas_chat(3, "lower split threshold: %.3lf keV\n", 
			xmldata->det->threshold_split_lo_keV);
	  }
	}

	else if (!strcmp(Uelement, "THRESHOLD_SPLIT_LO_FRACTION")) {
	  if (!strcmp(Uattribute, "VALUE")) {
	    xmldata->det->threshold_split_lo_fraction = (float)atof(attr[i+1]);
	    headas_chat(3, "lower split threshold: %.1lf %%\n", 
			xmldata->det->threshold_split_lo_fraction*100.);
	  }
	}
      
	if (EXIT_SUCCESS!=xmldata->status) return;
      } 
      // END of loop over different attributes.
    }
    // END of elements with independent attributes
  }
  // END of elements with attributes.
}



static void GenDetXMLElementEnd(void *data, const char *el) 
{
  struct XMLData* xmldata = (struct XMLData*)data;

  (void)el; // Unused parameter.

  // Check if an error has occurred previously.
  if (EXIT_SUCCESS!=xmldata->status) return;

  return;
}



void addGenDetPhotonImpact(GenDet* const det, const Impact* const impact, 
			   int* const status)
{

  // Call the detector operating clock routine.
  operateGenDetClock(det, impact->time, status);
  if (EXIT_SUCCESS!=*status) return;

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
    return; // Break the function (photon is not detected).
  }

  // Get the corresponding created charge.
  // NOTE: In this simulation the charge is represented by the nominal
  // photon energy [keV] which corresponds to the PHA channel according 
  // to the EBOUNDS table.
  float charge = getEnergy(channel, det->rmf, 0);
  assert(charge>=0.);

  // Create split events.
  makeGenSplitEvents(det->split, &impact->position, charge, det->pixgrid, det->line);
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
    getClockListElement(det->clocklist, time, &type, &element);
    switch (type) {
    case CL_NONE:
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: Clock list is empty!\n", *status);
      return;
    case CL_WAIT:
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
  } while(type!=CL_WAIT);
}



void GenDetLineShift(GenDet* const det)
{
  int ii;
  headas_chat(5, "lineshift\n");

  // Check if the detector contains more than 1 line.
  if (2>det->pixgrid->ywidth) return;

  // Apply the Charge Transfer Efficiency.
  if (det->cte!=1.) {
    int jj;
    for (ii=1; ii<det->pixgrid->ywidth; ii++) {
      if (0!=det->line[ii]->anycharge) {
	for (jj=0; jj<det->line[ii]->xwidth; jj++) {
	  det->line[ii]->charge[jj] *= det->cte;
	}
      }
    }
  }

  // Add the charges in line 1 to line 0.
  addGenDetLine(det->line[0], det->line[1]);

  // Clear the charges in line 1, as they are now contained in line 0.
  clearGenDetLine(det->line[1]);

  // Switch the other lines in increasing order such that the cleared 
  // original line number 1 will end up as the last line.
  for (ii=1; ii<det->pixgrid->ywidth-1; ii++) {
    switchGenDetLines(&det->line[ii], &det->line[ii+1]);
  }
}



void GenDetReadoutLine(GenDet* const det, const int lineindex, 
		       const int readoutindex, int* const status)
{
  headas_chat(5, "read out line %d as %d\n", lineindex, readoutindex);

  // Event data structure.
  GenEvent event = {.time = 0.};

  while (readoutGenDetLine(det->line[lineindex], &event)) {

    // Apply the charge thresholds.
    if (event.charge<=det->threshold_readout_lo_keV) continue;

    if (det->threshold_readout_up_keV >= 0.) {
      if (event.charge>=det->threshold_readout_up_keV) continue;
    }

    // Apply the detector response.
    event.pha = getChannel(event.charge, det->rmf);

    // Store the additional information.
    event.rawy  = readoutindex;
    event.time  = det->clocklist->time;  // Time of detection.
    event.frame = det->clocklist->frame; // Frame of detection.
    event.pat_type = 0;
    event.pat_id   = 0;
    event.pat_alig = 0;

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
    //      if ((status = PILGetFname("fits_templates", parameters->eventlist_template))) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not read environment variable 'SIXT_FITS_TEMPLATES'!\n", 
		   *status);
    return;
    //      }
  }

  // Append the filename of the template file itself.
  strcat(template, "/");
  strcat(template, det->eventfile_template);

  // Open a new event file from the specified template.
  det->eventfile = openNewGenEventFile(filename, template, status);
}


