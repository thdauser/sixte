#include "geninst.h"


////////////////////////////////////////////////////////////////////
// Local data type declarations.
////////////////////////////////////////////////////////////////////


/** Data structure given to the XML handler to transfer data. */
struct XMLParseData {
  GenInst* inst;
  int status;
};


////////////////////////////////////////////////////////////////////
// Static variables.
////////////////////////////////////////////////////////////////////


/** Flag indicating that the erodetbkgrndgen module is initialized and
    operational. */
static int eroBkgInitialized=0;


////////////////////////////////////////////////////////////////////
// Static function declarations.
////////////////////////////////////////////////////////////////////


/** Handler for the start of an XML element. */
static void GenInstXMLElementStart(void* data, const char* el, 
				   const char** attr);
/** Handler for the end of an XML element. */
static void GenInstXMLElementEnd(void* data, const char* el);


////////////////////////////////////////////////////////////////////
// Function Declarations.
////////////////////////////////////////////////////////////////////


void parseGenInstXML(GenInst* const inst, 
		     const char* const filename, 
		     int* const status);


////////////////////////////////////////////////////////////////////
// Program Code.
////////////////////////////////////////////////////////////////////


GenInst* newGenInst(int* const status) 
{
  // Allocate memory.
  GenInst* inst=(GenInst*)malloc(sizeof(GenInst));
  if (NULL==inst) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for GenInst failed");
    return(inst);
  }

  // Initialize all pointers with NULL.
  inst->tel=NULL;
  inst->det=NULL;
  inst->filename=NULL;
  inst->filepath=NULL;

  // Allocate memory for the GenTel and GenDet data structs.
  inst->det=newGenDet(status);
  CHECK_STATUS_RET(*status, inst);
  inst->tel=newGenTel(status);
  CHECK_STATUS_RET(*status, inst);

  return(inst);
}


void destroyGenInst(GenInst** const inst, int* const status)
{
  if (NULL!=*inst) {
    if (NULL!=(*inst)->tel) {
      destroyGenTel(&(*inst)->tel);
    }
    if (NULL!=(*inst)->det) {
      destroyGenDet(&(*inst)->det);
    }
    if (1==eroBkgInitialized) {
      eroBkgCleanUp(status);
      eroBkgInitialized=0;
    }
    if (NULL!=(*inst)->filename) {
      free((*inst)->filename);
    }
    if (NULL!=(*inst)->filepath) {
      free((*inst)->filepath);
    }
    free(*inst);
    *inst=NULL;
  }
}


GenInst* loadGenInst(const char* const filename, int* const status)
{
  // Get a new and empty data structure.
  GenInst* inst=newGenInst(status);
  CHECK_STATUS_RET(*status, inst);

  // Split the reference to the XML detector definition file
  // into path and filename. This has to be done before
  // calling the parser routine for the XML file.
  char filename2[MAXFILENAME];
  char rootname[MAXFILENAME];
  // Make a local copy of the filename variable in order to avoid
  // compiler warnings due to discarded const qualifier at the 
  // subsequent function call.
  strcpy(filename2, filename);
  fits_parse_rootname(filename2, rootname, status);
  CHECK_STATUS_RET(*status, inst);

  // Split rootname into the file path and the file name.
  char* lastslash=strrchr(rootname, '/');
  if (NULL==lastslash) {
    inst->filepath=(char*)malloc(sizeof(char));
    CHECK_NULL_RET(inst->filepath, *status, 
		   "memory allocation for filepath failed", inst);
    inst->filename=(char*)malloc((strlen(rootname)+1)*sizeof(char));
    CHECK_NULL_RET(inst->filename, *status, 
		   "memory allocation for filename failed", inst);
    strcpy(inst->filepath, "");
    strcpy(inst->filename, rootname);
  } else {
    lastslash++;
    inst->filename=(char*)malloc((strlen(lastslash)+1)*sizeof(char));
    CHECK_NULL_RET(inst->filename, *status, 
		   "memory allocation for filename failed", inst);
    strcpy(inst->filename, lastslash);
      
    *lastslash='\0';
    inst->filepath=(char*)malloc((strlen(rootname)+1)*sizeof(char));
    CHECK_NULL_RET(inst->filepath, *status, 
		   "memory allocation for filepath failed", inst);
    strcpy(inst->filepath, rootname);
  }
  // END of storing the filename and filepath.


  // Read in the XML definition of the detector.
  parseGenInstXML(inst, filename, status);
  CHECK_STATUS_RET(*status, inst);


  // Allocate memory for the detector pixels.
  inst->det->line=
    (GenDetLine**)malloc(inst->det->pixgrid->ywidth*sizeof(GenDetLine*));
  if (NULL==inst->det->line) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for GenDet pixel array failed");
    return(inst);
  }
  int ii;
  for (ii=0; ii<inst->det->pixgrid->ywidth; ii++) {
    inst->det->line[ii]=newGenDetLine(inst->det->pixgrid->xwidth, status);
    if (EXIT_SUCCESS!=*status) return(inst);
  }


  return(inst);
}


void parseGenInstXML(GenInst* const inst, 
		     const char* const filename, 
		     int* const status)
{
  headas_chat(5, "read instrument setup from XML file '%s' ...\n", filename);

  // Read the XML data from the file.
  // Open the specified file.
  FILE* xmlfile=fopen(filename, "r");
  if (NULL==xmlfile) {
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "failed opening XML "
	    "file '%s' for read access", filename);
    SIXT_ERROR(msg);
    return;
  }

  // The data is read from the XML file and stored in xmlbuffer
  // without any modifications.
  struct XMLBuffer* xmlbuffer=newXMLBuffer(status);
  CHECK_STATUS_VOID(*status);

  // Input buffer with an additional byte at the end for the 
  // termination of the string.
  const int buffer_size=256;
  char buffer[buffer_size+1];
  // Number of chars in buffer.
  int len;

  // Read all data from the file.
  do {
    // Get a piece of input into the buffer.
    len=fread(buffer, 1, buffer_size, xmlfile);
    buffer[len]='\0'; // Terminate the string.
    addString2XMLBuffer(xmlbuffer, buffer, status);
    CHECK_STATUS_VOID(*status);
  } while (!feof(xmlfile));

  // Close the file handler to the XML file.
  fclose(xmlfile);


  // Before acutally parsing the XML code, expand the loops and 
  // arithmetic operations in the GenDet XML description.
  // The expansion algorithm repeatetly scans the XML code and
  // searches for loop tags. It replaces the loop tags by repeating
  // the contained XML code.
  expandXML(xmlbuffer, status);
  CHECK_STATUS_VOID(*status);


  // Parse XML code in the xmlbuffer using the expat library.
  // Get an XML_Parser object.
  XML_Parser parser=XML_ParserCreate(NULL);
  if (NULL==parser) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("could not allocate memory for XML parser");
    return;
  }

  // Set data that is passed to the handler functions.
  struct XMLParseData xmlparsedata = {
    .inst   = inst,
    .status = EXIT_SUCCESS
  };
  XML_SetUserData(parser, &xmlparsedata);

  // Set the handler functions.
  XML_SetElementHandler(parser, GenInstXMLElementStart, GenInstXMLElementEnd);

  // Parse all the data in the string buffer.
  const int done=1;
  if (!XML_Parse(parser, xmlbuffer->text, strlen(xmlbuffer->text), done)) {
    // Parse error.
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "faild parsing XML file '%s':\n%s\n", 
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
  if (0==inst->det->pixgrid->xwidth) {
    headas_printf("*** warning: no specification of x-width of pixel array\n");
  }  
  if (0==inst->det->pixgrid->ywidth) {
    headas_printf("*** warning: no specification of y-width of pixel array\n");
  }

  if (0.==inst->det->pixgrid->xrpix) {
    headas_printf("*** warning: no specification of x reference pixel\n");
  }
  if (0.==inst->det->pixgrid->yrpix) {
    headas_printf("*** warning: no specification of y reference pixel\n");
  }

  if (0.==inst->det->pixgrid->xdelt) {
    headas_printf("*** warning: no specification of pixel x-width\n");
  }
  if (0.==inst->det->pixgrid->ydelt) {
    headas_printf("*** warning: no specification of pixel y-width\n");
  }
  
  if (0.>inst->det->pixgrid->xborder) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("invalid specification of x-border of pixels");
    return;    
  }
  if (0.>inst->det->pixgrid->yborder) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("invalid specification of y-border of pixels");
    return;    
  }

  if (NULL==inst->det->rmf) {
    headas_printf("*** warning: no specification of response file (RMF/RSP)\n");
  }
  if (NULL==inst->tel->arf) {
    headas_printf("*** warning: no specification of ARF\n");
  }

  if (0.==inst->tel->focal_length) {
    headas_printf("*** warning: no specification of the focal length of the telescope\n");
  }
  if (0.==inst->tel->fov_diameter) {
    headas_printf("*** warning: no specification of the diameter of the telescope FoV\n");
  }

  if (0==inst->det->readout_trigger) {
    headas_printf("*** warning: no specification of the readout trigger\n");
  }

  if (GS_EXPONENTIAL==inst->det->split->type) {
    if (inst->det->split->par1==0.) {
      *status=EXIT_FAILURE;
      SIXT_ERROR("no valid split model parameters in the XML file");
      return;    
    }
  }

  if (GS_GAUSS==inst->det->split->type) {
    if ((inst->det->split->par1==0.)&&(inst->det->split->par2==0.)) {
      *status=EXIT_FAILURE;
      SIXT_ERROR("no valid split model parameters in the XML file");
      return;    
    }
  }
  // END of checking, if all detector parameters have successfully been 
  // read from the XML file.
}


static void GenInstXMLElementStart(void* parsedata, 
				   const char* el, 
				   const char** attr) 
{
  struct XMLParseData* xmlparsedata=(struct XMLParseData*)parsedata;
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
    xmlparsedata->inst->tel->telescope=
      (char*)malloc((strlen(telescope)+1)*sizeof(char));
    CHECK_NULL_VOID(xmlparsedata->inst->tel->telescope, 
		      xmlparsedata->status,
		      "memory allocation for telescope name failed");
    strcpy(xmlparsedata->inst->tel->telescope, telescope);

  } else if (!strcmp(Uelement, "LINESHIFT")) {
    CLLineShift* cllineshift=newCLLineShift(&xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);
    append2ClockList(xmlparsedata->inst->det->clocklist, CL_LINESHIFT, 
		     cllineshift, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

  } else if (!strcmp(Uelement, "NEWFRAME")) {
    CLNewFrame* clnewframe=newCLNewFrame(&xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);
    append2ClockList(xmlparsedata->inst->det->clocklist, CL_NEWFRAME, 
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
    append2ClockList(xmlparsedata->inst->det->clocklist, CL_READOUTLINE, 
		     clreadoutline, &xmlparsedata->status);
      
  } else if (!strcmp(Uelement, "DIMENSIONS")) {

    xmlparsedata->inst->det->pixgrid->xwidth=
      getXMLAttributeInt(attr, "XWIDTH");
    xmlparsedata->inst->det->pixgrid->ywidth=
      getXMLAttributeInt(attr, "YWIDTH");
      
  } else if (!strcmp(Uelement, "WCS")) {

    xmlparsedata->inst->det->pixgrid->xrpix=
      getXMLAttributeFloat(attr, "XRPIX");
    xmlparsedata->inst->det->pixgrid->yrpix=
      getXMLAttributeFloat(attr, "YRPIX");
    xmlparsedata->inst->det->pixgrid->xrval=
      getXMLAttributeFloat(attr, "XRVAL");
    xmlparsedata->inst->det->pixgrid->yrval=
      getXMLAttributeFloat(attr, "YRVAL");
    xmlparsedata->inst->det->pixgrid->xdelt=
      getXMLAttributeFloat(attr, "XDELT");
    xmlparsedata->inst->det->pixgrid->ydelt=
      getXMLAttributeFloat(attr, "YDELT");
    xmlparsedata->inst->det->pixgrid->rota=
      getXMLAttributeFloat(attr, "ROTA")*M_PI/180.;
	
  } else if (!strcmp(Uelement, "PIXELBORDER")) {

    xmlparsedata->inst->det->pixgrid->xborder=getXMLAttributeFloat(attr, "X");
    xmlparsedata->inst->det->pixgrid->yborder=getXMLAttributeFloat(attr, "Y");
    
  } else if (!strcmp(Uelement, "RMF")) {

    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->inst->filepath);
    strcat(filepathname, filename);
    xmlparsedata->inst->det->rmf=loadRMF(filepathname, &xmlparsedata->status);

  } else if (!strcmp(Uelement, "ARF")) {

    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->inst->filepath);
    strcat(filepathname, filename);
    xmlparsedata->inst->tel->arf=loadARF(filepathname, &xmlparsedata->status);

  } else if (!strcmp(Uelement, "PSF")) {
    
    // The focal length must be specified before load the PSF.
    // Check if this is the case.
    if (xmlparsedata->inst->tel->focal_length<=0.) {
      xmlparsedata->status=EXIT_FAILURE;
      SIXT_ERROR("telescope focal length must be specified "
		 "before loading the PSF");
      return;
    }
    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->inst->filepath);
    strcat(filepathname, filename);
    xmlparsedata->inst->tel->psf= 
      newPSF(filepathname, xmlparsedata->inst->tel->focal_length, 
	     &xmlparsedata->status);

  } else if (!strcmp(Uelement, "CODEDMASK")) {

    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->inst->filepath);
    strcat(filepathname, filename);
    xmlparsedata->inst->tel->coded_mask = 
      getCodedMaskFromFile(filepathname, &xmlparsedata->status);

  } else if (!strcmp(Uelement, "VIGNETTING")) {

    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->inst->filepath);
    strcat(filepathname, filename);
    xmlparsedata->inst->tel->vignetting =
      newVignetting(filepathname, &xmlparsedata->status);

  } else if (!strcmp(Uelement, "FOCALLENGTH")) {

    xmlparsedata->inst->tel->focal_length=getXMLAttributeFloat(attr, "VALUE");

  } else if (!strcmp(Uelement, "FOV")) {

    xmlparsedata->inst->tel->fov_diameter=
      getXMLAttributeFloat(attr, "DIAMETER")*M_PI/180.;

  } else if (!strcmp(Uelement, "CTE")) {

    xmlparsedata->inst->det->cte=getXMLAttributeFloat(attr, "VALUE");
	
  } else if (!strcmp(Uelement, "BADPIXMAP")) {

    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->inst->filepath);
    strcat(filepathname, filename);
    xmlparsedata->inst->det->badpixmap = 
      loadBadPixMap(filepathname, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

  } else if (!strcmp(Uelement, "EROBACKGROUND")) {

    if (0==eroBkgInitialized) {
      char filename[MAXFILENAME];
      getXMLAttributeString(attr, "FILENAME", filename);
      char filepathname[MAXFILENAME];
      strcpy(filepathname, xmlparsedata->inst->filepath);
      strcat(filepathname, filename);
      eroBkgInitialize(filepathname, &xmlparsedata->status);
      CHECK_STATUS_VOID(xmlparsedata->status);
      eroBkgInitialized=1;
    }

    xmlparsedata->inst->det->erobackground=1;

  } else if (!strcmp(Uelement, "SPLIT")) {

    char type[MAXMSG];
    getXMLAttributeString(attr, "TYPE", type);
    strtoupper(type);
    if (!strcmp(type, "NONE")) {
      xmlparsedata->inst->det->split->type = GS_NONE;
    } else if (!strcmp(type, "GAUSS")) {
      xmlparsedata->inst->det->split->type = GS_GAUSS;
    } else if (!strcmp(type, "EXPONENTIAL")) {
      xmlparsedata->inst->det->split->type = GS_EXPONENTIAL;
    }
    xmlparsedata->inst->det->split->par1=getXMLAttributeFloat(attr, "PAR1");
    xmlparsedata->inst->det->split->par2=getXMLAttributeFloat(attr, "PAR2");

  } else if (!strcmp(Uelement, "READOUT")) {

    char mode[MAXMSG];
    getXMLAttributeString(attr, "MODE", mode);
    strtoupper(mode);
    if (!strcmp(mode, "TIME")) {
      xmlparsedata->inst->det->readout_trigger=GENDET_TIME_TRIGGERED;
    } else if (!strcmp(mode, "EVENT")) {
      xmlparsedata->inst->det->readout_trigger=GENDET_EVENT_TRIGGERED;
    }
      
  } else if (!strcmp(Uelement, "WAIT")) {
    
    float waittime=getXMLAttributeFloat(attr, "TIME");
    CLWait* clwait=newCLWait(waittime, &xmlparsedata->status);
    append2ClockList(xmlparsedata->inst->det->clocklist, CL_WAIT, 
		     clwait, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

    // Accumulate the amount of time required for one read-out frame.
    xmlparsedata->inst->det->frametime += waittime;
	
  } else if (!strcmp(Uelement, "CLEARLINE")) {

    int lineindex=getXMLAttributeInt(attr, "LINEINDEX");
    CLClearLine* clclearline=newCLClearLine(lineindex, &xmlparsedata->status);
    append2ClockList(xmlparsedata->inst->det->clocklist, CL_CLEARLINE, 
		     clclearline, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

  } else if (!strcmp(Uelement, "THRESHOLD_READOUT_LO_KEV")) {

    xmlparsedata->inst->det->threshold_readout_lo_keV = 
      getXMLAttributeFloat(attr, "VALUE");
    headas_chat(3, "lower readout threshold: %.3lf keV\n", 
		xmlparsedata->inst->det->threshold_readout_lo_keV);
	
  } else if (!strcmp(Uelement, "THRESHOLD_PATTERN_UP_KEV")) {

    xmlparsedata->inst->det->threshold_pattern_up_keV = 
      getXMLAttributeFloat(attr, "VALUE");
    headas_chat(3, "upper pattern threshold: %.3lf keV\n", 
		xmlparsedata->inst->det->threshold_pattern_up_keV);

  } else if (!strcmp(Uelement, "THRESHOLD_EVENT_LO_KEV")) {

    xmlparsedata->inst->det->threshold_event_lo_keV = 
      getXMLAttributeFloat(attr, "VALUE");
    headas_chat(3, "lower event threshold: %.3lf keV\n", 
		xmlparsedata->inst->det->threshold_event_lo_keV);

  } else if (!strcmp(Uelement, "THRESHOLD_SPLIT_LO_KEV")) {

    xmlparsedata->inst->det->threshold_split_lo_keV =
      getXMLAttributeFloat(attr, "VALUE");
    headas_chat(3, "lower split threshold: %.3lf keV\n", 
		xmlparsedata->inst->det->threshold_split_lo_keV);

  } else if (!strcmp(Uelement, "THRESHOLD_SPLIT_LO_FRACTION")) {

    xmlparsedata->inst->det->threshold_split_lo_fraction = 
      getXMLAttributeFloat(attr, "VALUE");
    headas_chat(3, "lower split threshold: %.1lf %%\n", 
		xmlparsedata->inst->det->threshold_split_lo_fraction*100.);

  }
  CHECK_STATUS_VOID(xmlparsedata->status);
}


static void GenInstXMLElementEnd(void* parsedata, const char* el) 
{
  struct XMLParseData* xmlparsedata=(struct XMLParseData*)parsedata;

  (void)el; // Unused parameter.

  // Check if an error has occurred previously.
  CHECK_STATUS_VOID(xmlparsedata->status);

  return;
}

