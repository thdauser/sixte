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
  det->pixgrid->xwidth = 0;
  det->pixgrid->ywidth = 0;
  det->pixgrid->xrpix  = 0.;
  det->pixgrid->yrpix  = 0.;
  det->pixgrid->xrval  = 0.;
  det->pixgrid->yrval  = 0.;
  det->pixgrid->xdelt  = 0.;
  det->pixgrid->ydelt  = 0.;
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
    SIXT_ERROR("could not allocate memory for XML parser");
    return;
  }

  // Set data that is passed to the handler functions.
  struct XMLParseData xmlparsedata = {
    .det    = det,
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
    SIXT_ERROR(msg);
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
  if (0==det->pixgrid->xwidth) {
    headas_printf("*** warning: no specification of x-width of pixel array\n");
  }  
  if (0==det->pixgrid->ywidth) {
    headas_printf("*** warning: no specification of y-width of pixel array\n");
  }

  if (0.==det->pixgrid->xrpix) {
    headas_printf("*** warning: no specification of x reference pixel\n");
  }
  if (0.==det->pixgrid->yrpix) {
    headas_printf("*** warning: no specification of y reference pixel\n");
  }

  if (0.==det->pixgrid->xrval) {
    headas_printf("*** warning: no specification of x reference value\n");
  }
  if (0.==det->pixgrid->yrval) {
    headas_printf("*** warning: no specification of y reference value\n");
  }

  if (0.==det->pixgrid->xdelt) {
    headas_printf("*** warning: no specification of pixel x-width\n");
  }
  if (0.==det->pixgrid->ydelt) {
    headas_printf("*** warning: no specification of pixel y-width\n");
  }
  
  if (0.>det->pixgrid->xborder) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("invalid specification of x-border of pixels");
    return;    
  }
  if (0.>det->pixgrid->yborder) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("invalid specification of y-border of pixels");
    return;    
  }

  if (NULL==det->rmf) {
    headas_printf("*** warning: no specification of response file (RMF/RSP)\n");
  }
  if (NULL==det->arf) {
    headas_printf("*** warning: no specification of ARF\n");
  }

  if (0.==det->focal_length) {
    headas_printf("*** warning: no specification of the focal length of the telescope\n");
  }
  if (0.==det->fov_diameter) {
    headas_printf("*** warning: no specification of the diameter of the telescope FoV\n");
  }

  if (0==det->readout_trigger) {
    headas_printf("*** warning: no specification of the readout trigger\n");
  }

  if (GS_EXPONENTIAL==det->split->type) {
    if (det->split->par1==0.) {
      *status = EXIT_FAILURE;
      SIXT_ERROR("no valid split model parameters in the XML file");
      return;    
    }
  }

  if (GS_GAUSS==det->split->type) {
    if ((det->split->par1==0.)&&(det->split->par2==0.)) {
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


static void GenDetXMLElementStart(void* parsedata, const char* el, 
				  const char** attr) 
{
  struct XMLParseData* xmlparsedata = (struct XMLParseData*)parsedata;
  char Uelement[MAXMSG]; // Upper case version of XML element.

  // Check if an error has occurred previously.
  CHECK_STATUS_VOID(xmlparsedata->status);

  // Convert the element to an upper case string.
  strcpy(Uelement, el);
  strtoupper(Uelement);

  // Check for different elements.
  if (!strcmp(Uelement, "TELESCOP")) {
    // Determine the name of the telescope.
    char telescope[MAXMSG];
    getXMLAttributeString(attr, "NAME", telescope);
    xmlparsedata->det->telescope=
      (char*)malloc((strlen(telescope)+1)*sizeof(char));
    CHECK_NULL_VOID(xmlparsedata->det->telescope, 
		      xmlparsedata->status,
		      "memory allocation for telescope name failed");
    strcpy(xmlparsedata->det->telescope, telescope);

  } else if (!strcmp(Uelement, "LINESHIFT")) {
    CLLineShift* cllineshift=newCLLineShift(&xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);
    append2ClockList(xmlparsedata->det->clocklist, CL_LINESHIFT, 
		     cllineshift, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

  } else if (!strcmp(Uelement, "NEWFRAME")) {
    CLNewFrame* clnewframe=newCLNewFrame(&xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);
    append2ClockList(xmlparsedata->det->clocklist, CL_NEWFRAME, 
		     clnewframe, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

  } else if (!strcmp(Uelement, "READOUTLINE")) {

    int lineindex=getXMLAttributeInt(attr, "LINEINDEX");
    if (lineindex<0) {
      xmlparsedata->status=EXIT_FAILURE;
      SIXT_ERROR("negative index for readout line");
      return;
    }
    int readoutindex=getXMLAttributeInt(attr, "READOUTINDEX");
    if (readoutindex<0) {
      xmlparsedata->status=EXIT_FAILURE;
      SIXT_ERROR("negative index for readout line");
      return;
    }
    CLReadoutLine* clreadoutline=
      newCLReadoutLine(lineindex, readoutindex, &xmlparsedata->status);
    append2ClockList(xmlparsedata->det->clocklist, CL_READOUTLINE, 
		     clreadoutline, &xmlparsedata->status);
      
  } else if (!strcmp(Uelement, "DIMENSIONS")) {

    xmlparsedata->det->pixgrid->xwidth=getXMLAttributeInt(attr, "XWIDTH");
    xmlparsedata->det->pixgrid->ywidth=getXMLAttributeInt(attr, "YWIDTH");
      
  } else if (!strcmp(Uelement, "WCS")) {

    xmlparsedata->det->pixgrid->xrpix=getXMLAttributeFloat(attr, "XRPIX");
    xmlparsedata->det->pixgrid->yrpix=getXMLAttributeFloat(attr, "YRPIX");
    xmlparsedata->det->pixgrid->xrval=getXMLAttributeFloat(attr, "XRVAL");
    xmlparsedata->det->pixgrid->yrval=getXMLAttributeFloat(attr, "YRVAL");
    xmlparsedata->det->pixgrid->xdelt=getXMLAttributeFloat(attr, "XDELT");
    xmlparsedata->det->pixgrid->ydelt=getXMLAttributeFloat(attr, "YDELT");
	
  } else if (!strcmp(Uelement, "PIXELBORDER")) {

    xmlparsedata->det->pixgrid->xborder=getXMLAttributeFloat(attr, "X");
    xmlparsedata->det->pixgrid->yborder=getXMLAttributeFloat(attr, "Y");
    
  } else if (!strcmp(Uelement, "RMF")) {

    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->det->filepath);
    strcat(filepathname, filename);
    xmlparsedata->det->rmf=loadRMF(filepathname, &xmlparsedata->status);

  } else if (!strcmp(Uelement, "ARF")) {

    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->det->filepath);
    strcat(filepathname, filename);
    xmlparsedata->det->arf=loadARF(filepathname, &xmlparsedata->status);

  } else if (!strcmp(Uelement, "PSF")) {
    
    // The focal length must be specified before load the PSF.
    // Check if this is the case.
    if (xmlparsedata->det->focal_length<=0.) {
      xmlparsedata->status=EXIT_FAILURE;
      SIXT_ERROR("telescope focal length must be specified "
		 "before loading the PSF");
      return;
    }
    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->det->filepath);
    strcat(filepathname, filename);
    xmlparsedata->det->psf= 
      newPSF(filepathname, xmlparsedata->det->focal_length, 
	     &xmlparsedata->status);

  } else if (!strcmp(Uelement, "CODEDMASK")) {

    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->det->filepath);
    strcat(filepathname, filename);
    xmlparsedata->det->coded_mask = 
      getCodedMaskFromFile(filepathname, &xmlparsedata->status);

  } else if (!strcmp(Uelement, "VIGNETTING")) {

    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->det->filepath);
    strcat(filepathname, filename);
    xmlparsedata->det->vignetting =
      newVignetting(filepathname, &xmlparsedata->status);

  } else if (!strcmp(Uelement, "FOCALLENGTH")) {
    
    xmlparsedata->det->focal_length=getXMLAttributeFloat(attr, "VALUE");

  } else if (!strcmp(Uelement, "FOV")) {

    xmlparsedata->det->fov_diameter=
      getXMLAttributeFloat(attr, "DIAMETER")*M_PI/180.;

  } else if (!strcmp(Uelement, "CTE")) {
    xmlparsedata->det->cte=getXMLAttributeFloat(attr, "VALUE");
	
  } else if (!strcmp(Uelement, "BADPIXMAP")) {

    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->det->filepath);
    strcat(filepathname, filename);
    xmlparsedata->det->badpixmap = 
      loadBadPixMap(filepathname, &xmlparsedata->status);

  } else if (!strcmp(Uelement, "EROBACKGROUND")) {

    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->det->filepath);
    strcat(filepathname, filename);
    eroBkgInitialize(filepathname, &xmlparsedata->status);
    xmlparsedata->det->erobackground=1;

  } else if (!strcmp(Uelement, "SPLIT")) {

    char type[MAXMSG];
    getXMLAttributeString(attr, "TYPE", type);
    strtoupper(type);
    if (!strcmp(type, "NONE")) {
      xmlparsedata->det->split->type = GS_NONE;
    } else if (!strcmp(type, "GAUSS")) {
      xmlparsedata->det->split->type = GS_GAUSS;
    } else if (!strcmp(type, "EXPONENTIAL")) {
      xmlparsedata->det->split->type = GS_EXPONENTIAL;
    }
    xmlparsedata->det->split->par1=getXMLAttributeFloat(attr, "PAR1");
    xmlparsedata->det->split->par2=getXMLAttributeFloat(attr, "PAR2");

  } else if (!strcmp(Uelement, "READOUT")) {

    char mode[MAXMSG];
    getXMLAttributeString(attr, "MODE", mode);
    strtoupper(mode);
    if (!strcmp(mode, "TIME")) {
      xmlparsedata->det->readout_trigger=GENDET_TIME_TRIGGERED;
    } else if (!strcmp(mode, "EVENT")) {
      xmlparsedata->det->readout_trigger=GENDET_EVENT_TRIGGERED;
    }
      
  } else if (!strcmp(Uelement, "WAIT")) {
    
    float waittime=getXMLAttributeFloat(attr, "TIME");
    CLWait* clwait=newCLWait(waittime, &xmlparsedata->status);
    append2ClockList(xmlparsedata->det->clocklist, CL_WAIT, 
		     clwait, &xmlparsedata->status);

    // Accumulate the amount of time required for one read-out frame.
    xmlparsedata->det->frametime += waittime;
	
  } else if (!strcmp(Uelement, "CLEARLINE")) {

    int lineindex=getXMLAttributeInt(attr, "LINEINDEX");
    CLClearLine* clclearline=newCLClearLine(lineindex, &xmlparsedata->status);
    append2ClockList(xmlparsedata->det->clocklist, CL_CLEARLINE, 
		     clclearline, &xmlparsedata->status);

  } else if (!strcmp(Uelement, "THRESHOLD_READOUT_LO_KEV")) {

    xmlparsedata->det->threshold_readout_lo_keV = 
      getXMLAttributeFloat(attr, "VALUE");
    headas_chat(3, "lower readout threshold: %.3lf keV\n", 
		xmlparsedata->det->threshold_readout_lo_keV);
	
  } else if (!strcmp(Uelement, "THRESHOLD_READOUT_UP_KEV")) {

    xmlparsedata->det->threshold_readout_up_keV = 
      getXMLAttributeFloat(attr, "VALUE");
    headas_chat(3, "upper readout threshold: %.3lf keV\n", 
		xmlparsedata->det->threshold_readout_up_keV);

  } else if (!strcmp(Uelement, "THRESHOLD_READOUT_LO_PHA")) {

    xmlparsedata->det->threshold_readout_lo_PHA = 
      getXMLAttributeLong(attr, "VALUE");
    headas_chat(3, "lower readout threshold: %ld PHA\n", 
		xmlparsedata->det->threshold_readout_lo_PHA);

  } else if (!strcmp(Uelement, "THRESHOLD_READOUT_UP_PHA")) {

    xmlparsedata->det->threshold_readout_up_PHA = 
      getXMLAttributeLong(attr, "VALUE");
    headas_chat(3, "upper readout threshold: %ld PHA\n", 
		xmlparsedata->det->threshold_readout_up_PHA);

  } else if (!strcmp(Uelement, "THRESHOLD_EVENT_LO_KEV")) {

    xmlparsedata->det->threshold_event_lo_keV = 
      getXMLAttributeFloat(attr, "VALUE");
    headas_chat(3, "lower event threshold: %.3lf keV\n", 
		xmlparsedata->det->threshold_event_lo_keV);

  } else if (!strcmp(Uelement, "THRESHOLD_SPLIT_LO_KEV")) {

    xmlparsedata->det->threshold_split_lo_keV =
      getXMLAttributeFloat(attr, "VALUE");
    headas_chat(3, "lower split threshold: %.3lf keV\n", 
		xmlparsedata->det->threshold_split_lo_keV);

  } else if (!strcmp(Uelement, "THRESHOLD_SPLIT_LO_FRACTION")) {

    xmlparsedata->det->threshold_split_lo_fraction = 
      getXMLAttributeFloat(attr, "VALUE");
    headas_chat(3, "lower split threshold: %.1lf %%\n", 
		xmlparsedata->det->threshold_split_lo_fraction*100.);

  }
  if (EXIT_SUCCESS!=xmlparsedata->status) return;
}


static void GenDetXMLElementEnd(void* parsedata, const char* el) 
{
  struct XMLParseData* xmlparsedata = (struct XMLParseData*)parsedata;

  (void)el; // Unused parameter.

  // Check if an error has occurred previously.
  if (EXIT_SUCCESS!=xmlparsedata->status) return;

  return;
}


