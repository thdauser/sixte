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


   Copyright 2007 - 2018: Christian Schmid, Mirjam Oertel, FAU.
   Manuel Castro, National Institute for Space Research (INPE),
		 Brazil; under grant #2017/00968-6,
		 São Paulo Research Foundation (FAPESP).
     
   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg

*/

#include "masksystem.h"


////////////////////////////////////////////////////////////////////
// Local data type declarations.
////////////////////////////////////////////////////////////////////


/** Data structure given to the XML handler to transfer data. */
struct XMLParseData {
  MaskSystem* inst;
  unsigned int seed;
  int status;
};


////////////////////////////////////////////////////////////////////
// Static variables.
////////////////////////////////////////////////////////////////////


/** Flag indicating that the erodetbkgrndgen module is initialized and
    operational. */
static int auxBkgInitialized=0;


////////////////////////////////////////////////////////////////////
// Static function declarations.
////////////////////////////////////////////////////////////////////


/** Handler for the start of an XML element. */
static void MaskSystemXMLElementStart(void* data, const char* el,
				   const char** attr);
/** Handler for the end of an XML element. */
static void MaskSystemXMLElementEnd(void* data, const char* el);



////////////////////////////////////////////////////////////////////
// Program Code.
////////////////////////////////////////////////////////////////////


MaskSystem* newMaskSystem(int* const status)
{
  // Allocate memory.
  MaskSystem* inst=(MaskSystem*)malloc(sizeof(MaskSystem));
  if (NULL==inst) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for MaskSystem failed");
    return(inst);
  }

  // Initialize all pointers with NULL.

  inst->filename=NULL;
  inst->filepath=NULL;
  inst->x_mask=0;
  inst->y_mask=0;
  inst->mask_distance=0;
  inst->mask_pattern=NULL;
  inst->collimator_height=0;
  inst->wall_thickness=0;
  inst->flag=0;
  inst->x_det=0;
  inst->y_det=0;
  inst->det_pixelwidth;
  inst->det_width=0;
  inst->width=0;
  inst->DCU_length=0;
  inst->DCU_gap=0;
  inst->DCA_gap=0;
  inst->pixelwidth=0;
  inst->repixsize=0;
 // // Allocate memory for the GenTel and GenDet data structs.
 // inst->det=newGenDet(status);
 // CHECK_STATUS_RET(*status, inst);
 // inst->tel=newGenTel(status);
 // CHECK_STATUS_RET(*status, inst);

  return(inst);
}


void destroyMaskSystem(MaskSystem** const inst, int* const status)
{
  if (NULL!=*inst) {
//    if (NULL!=(*inst)->tel) {
//      destroyGenTel(&(*inst)->tel);
//    }
//    if (NULL!=(*inst)->det) {
//      destroyGenDet(&(*inst)->det);
//    }
    if (NULL!=(*inst)->filename) {
      free((*inst)->filename);
    }
    if (NULL!=(*inst)->filepath) {
      free((*inst)->filepath);
    }
    if (1==auxBkgInitialized) {
      bkgCleanUp(status);
      auxBkgInitialized=0;
    }
    (*inst)->x_mask=0;
    (*inst)->y_mask=0;
    (*inst)->mask_distance=0;
    if (NULL!=(*inst)->mask_pattern) {
      free((*inst)->mask_pattern);
    }
    (*inst)->collimator_height=0;
    (*inst)->wall_thickness=0;
    (*inst)->flag=0;
    (*inst)->x_det=0;
    (*inst)->y_det=0;
    (*inst)->det_pixelwidth=0;
    (*inst)->det_width=0;
    (*inst)->width=0;
    (*inst)->DCU_length=0;
    (*inst)->DCU_gap=0;
    (*inst)->DCA_gap=0;
    (*inst)->pixelwidth=0;
    (*inst)->repixsize=0;


//    if (0!=(*inst)->det_pixelwidth) {
//      free((*inst)->det_pixelwidth);
//    }
//    if (0!=(*inst)->det_width) {
//      free((*inst)->det_width);
//    }
//    if (0!=(*inst)->width) {
//      free((*inst)->width);
//    }
//    if (0!=(*inst)->DCU_length) {
//      free((*inst)->DCU_length);
//    }
//    if (0!=(*inst)->DCU_gap) {
//      free((*inst)->DCU_gap);
//    }
//    if (0!=(*inst)->DCA_gap) {
//      free((*inst)->DCA_gap);
//    }
//    if (0!=(*inst)->pixelwidth) {
//      free((*inst)->pixelwidth);
//    }
//    if (0!=(*inst)->repixsize) {
//      free((*inst)->repixsize);
//    }
    free(*inst);
    *inst=NULL;
  }
}


void parseMaskSystemXML(MaskSystem* const inst,
		     const char* const filename,
		     const unsigned int seed,
		     int* const status)
{
  headas_chat(5, "read instrument setup from advanced XML file '%s' ...\n", filename);

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

  // The data are read from the XML file and stored in xmlbuffer
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

  // Before expanding loops in the XML file, add the included code to it.
  expandIncludesXML(xmlbuffer, filename, status);
  CHECK_STATUS_VOID(*status);

  // Before actually parsing the XML code, expand the loops and
  // arithmetic operations in the GenDet XML description.
  // The expansion algorithm repeatedly scans the XML code and
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
  struct XMLParseData xmlparsedata={
    .inst  =inst,
    .seed  =seed,
    .status=EXIT_SUCCESS
  };
  XML_SetUserData(parser, &xmlparsedata);

  // Set the handler functions.
  XML_SetElementHandler(parser, MaskSystemXMLElementStart, MaskSystemXMLElementEnd);

  // Parse all the data in the string buffer.
  const int done=1;
  if (!XML_Parse(parser, xmlbuffer->text, strlen(xmlbuffer->text), done)) {
    // Parse error.
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "failed parsing XML file '%s':\n%s\n",
	    filename, XML_ErrorString(XML_GetErrorCode(parser)));
    printf("%s", xmlbuffer->text);
    SIXT_ERROR(msg);
    return;
  }
  // Check for errors.
  if (EXIT_SUCCESS!=xmlparsedata.status) {
    *status=xmlparsedata.status;
    return;
  }

  // Release memory.
  XML_ParserFree(parser);

  // Remove the XML string buffer.
  freeXMLBuffer(&xmlbuffer);

//  // Check if all required parameters have been read successfully from
//  // the XML file.
//  if (INT_MAX==inst->det->pixgrid->xwidth) {
//    *status=EXIT_FAILURE;
//    SIXT_ERROR("no specification of x-width of pixel array");
//    return;
//  }
//  if (INT_MAX==inst->det->pixgrid->ywidth) {
//    *status=EXIT_FAILURE;
//    SIXT_ERROR("no specification of y-width of pixel array");
//    return;
//  }
//
//  if (isnan(inst->det->pixgrid->xrpix)) {
//    *status=EXIT_FAILURE;
//    SIXT_ERROR("no specification of x reference pixel");
//    return;
//  }
//  if (isnan(inst->det->pixgrid->yrpix)) {
//    *status=EXIT_FAILURE;
//    SIXT_ERROR("no specification of y reference pixel");
//    return;
//  }
//
//  if (isnan(inst->det->pixgrid->xdelt)) {
//    *status=EXIT_FAILURE;
//    SIXT_WARNING("no specification of pixel x-width");
//    return;
//  }
//  if (isnan(inst->det->pixgrid->ydelt)) {
//    *status=EXIT_FAILURE;
//    SIXT_WARNING("no specification of pixel y-width");
//    return;
//  }
//
//  if (inst->det->pixgrid->xborder<0.) {
//    *status=EXIT_FAILURE;
//    SIXT_ERROR("invalid specification of x-border of pixels");
//    return;
//  }
//  if (inst->det->pixgrid->yborder<0.) {
//    *status=EXIT_FAILURE;
//    SIXT_ERROR("invalid specification of y-border of pixels");
//    return;
//  }
//
//  if (NULL==inst->det->rmf) {
//    SIXT_WARNING("no specification of response file (RMF/RSP)");
//  }
//  if (NULL==inst->tel->arf) {
//    SIXT_WARNING("no specification of ARF");
//  }
//
//  if (NULL==inst->tel->psf) {
//    SIXT_WARNING("no specification of PSF");
//  }
//
//  if (0.==inst->tel->focal_length) {
//    SIXT_WARNING("no specification of the focal length of the telescope");
//  }
//  if (0.==inst->tel->fov_diameter) {
//    SIXT_WARNING("no specification of the diameter of the telescope FoV");
//  }
//
//  if (0==inst->det->readout_trigger) {
//    SIXT_WARNING("no specification of the readout trigger");
//  }
//
//  if (GS_EXPONENTIAL==inst->det->split->type) {
//    if (inst->det->split->par1==0.) {
//      *status=EXIT_FAILURE;
//      SIXT_ERROR("no valid split model parameters in the XML file");
//      return;
//    }
//  }
//
//  if (GS_GAUSS==inst->det->split->type) {
//    if ((inst->det->split->par1==0.)&&(inst->det->split->par2==0.)) {
//      *status=EXIT_FAILURE;
//      SIXT_ERROR("no valid split model parameters in the XML file");
//      return;
//    }
//  }
//
//
//
//  // change borders for event driven detectors
//
//  if(GENDET_TIME_TRIGGERED!=inst->det->readout_trigger){
//    inst->det->rawymin=0;
//    inst->det->rawymax=inst->det->pixgrid->ywidth-1;
//  }

  // END of checking, if all detector parameters have successfully been
  // read from the XML file.
}



MaskSystem* loadMaskSystem(const char* const filename,
		     const unsigned int seed,
		     int* const status)
{
  // Get a new and empty data structure.
  MaskSystem* inst=newMaskSystem(status);
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
  parseMaskSystemXML(inst, filename, seed, status);
  CHECK_STATUS_RET(*status, inst);

//  // Allocate memory for the detector pixels.
//  inst->det->line=
//    (GenDetLine**)malloc(inst->det->pixgrid->ywidth*sizeof(GenDetLine*));
//  if (NULL==inst->det->line) {
//    *status=EXIT_FAILURE;
//    SIXT_ERROR("memory allocation for GenDet pixel array failed");
//    return(inst);
//  }
//  int ii;
//  for (ii=0; ii<inst->det->pixgrid->ywidth; ii++) {
//    inst->det->line[ii]=newGenDetLine(inst->det->pixgrid->xwidth, status);
//    if (EXIT_SUCCESS!=*status) return(inst);
//  }


  return(inst);
}


static void MaskSystemXMLElementStart(void* parsedata,
				   const char* el,
				   const char** attr)
{
  struct XMLParseData* xmlparsedata=(struct XMLParseData*)parsedata;

  // Check if an error has occurred previously.
  CHECK_STATUS_VOID(xmlparsedata->status);

  // Convert the element to an upper case string.
  char Uelement[MAXMSG];
  strcpy(Uelement, el);
  strtoupper(Uelement);

 //Read detector parameters

  if (!strcmp(Uelement, "DETECTOR")) {

    xmlparsedata->inst->x_det = getXMLAttributeFloat(attr, "X_DET");
    xmlparsedata->inst->y_det = getXMLAttributeFloat(attr, "Y_DET");
    xmlparsedata->inst->det_pixelwidth = getXMLAttributeFloat(attr, "DET_PIXELWIDTH");
    xmlparsedata->inst->det_width = getXMLAttributeFloat(attr, "DET_WIDTH");
    xmlparsedata->inst->width = getXMLAttributeInt(attr, "WIDTH");
    xmlparsedata->inst->DCU_length = getXMLAttributeFloat(attr, "DCU_LENGTH");
    xmlparsedata->inst->DCU_gap = getXMLAttributeFloat(attr, "DCU_GAP");
    xmlparsedata->inst->DCA_gap = getXMLAttributeFloat(attr, "DCA_GAP");
    xmlparsedata->inst->pixelwidth = getXMLAttributeFloat(attr, "PIXELWIDTH");
    xmlparsedata->inst->repixsize = getXMLAttributeFloat(attr, "REPIXSIZE");

}

 //Read mask parameters

  if (!strcmp(Uelement, "CODED_MASK")) {

    xmlparsedata->inst->x_mask = getXMLAttributeFloat(attr, "X_MASK");
    xmlparsedata->inst->y_mask = getXMLAttributeFloat(attr, "Y_MASK");
    xmlparsedata->inst->mask_distance = getXMLAttributeFloat(attr, "MASK_DISTANCE");
// Mask pattern
    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "MASK_PATTERN", filename);

    // Store the file name of the ARF.
    xmlparsedata->inst->mask_pattern =  (char*)malloc((strlen(filename)+1)*sizeof(char));
    CHECK_NULL_VOID(xmlparsedata->inst->mask_pattern,
		    xmlparsedata->status,
		    "memory allocation for mask pattern name failed");
    strcpy(xmlparsedata->inst->mask_pattern, filename);


    xmlparsedata->inst->collimator_height = getXMLAttributeFloat(attr, "COLLIMATOR_HEIGHT");
    xmlparsedata->inst->wall_thickness = getXMLAttributeFloat(attr, "WALL_THICKNESS");
    xmlparsedata->inst->flag = getXMLAttributeInt(attr, "FLAG");

}

}


static void MaskSystemXMLElementEnd(void* parsedata, const char* el)
{
  struct XMLParseData* xmlparsedata=(struct XMLParseData*)parsedata;

  (void)el;

  // Check if an error has occurred previously.
  CHECK_STATUS_VOID(xmlparsedata->status);
}
