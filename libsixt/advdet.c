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


   Copyright 2014 Thorsten Brand, FAU
*/

#include "advdet.h"

/** Data structure given to the XML handler to transfer data. */
struct XMLParseData {
  AdvDet* det;
  int status;
};

////////////////////////////////////////////////////////////////////
// Static function declarations.
////////////////////////////////////////////////////////////////////

/** Handler for the start of an XML element. */
static void AdvDetXMLElementStart(void* parsedata, 
				   const char* el, 
				   const char** attr);
/** Handler for the end of an XML element. */
static void AdvDetXMLElementEnd(void* parsedata, const char* el);

////////////////////////////////////////////////////////////////////
// Program Code
////////////////////////////////////////////////////////////////////

AdvDet* newAdvDet(int* const status){
  
  // Allocate memory.
  AdvDet* det=(AdvDet*)malloc(sizeof(AdvDet));
  if (NULL==det) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for AdvDet failed");
    return(det);
  }

  // Initialize all pointers with NULL.
  det->pix=NULL;
  det->filename=NULL;
  det->filepath=NULL;
  det->sx=0.;
  det->sy=0.;
  det->npix=0;
  det->cpix=0;

  return(det);
}

void destroyAdvDet(AdvDet **det, int* const status){
  
  if(NULL!=(*det)){
    if(NULL!=(*det)->pix){
      free((*det)->pix);
    }
  }
}

int CheckAdvPixImpact(AdvPix pix, Impact *imp){
  
  // Calculate impact coordinates in respect to the 
  // pixel coordinate system
  double u, v;
  
  u = imp->position.x - pix.sx;
  v = imp->position.y - pix.sy;
  
  // Calculate half width and height of the rectangular pixel
  double deltu, deltv;
  deltu=pix.width/2.;
  deltv=pix.height/2.;
 
  // Check if the impact lies in the rectangular pixel.
  // Return 1 if yes, 0 if not.
  if((u >= -deltu) && (u <= deltu) && (v >= -deltv) && (v <= deltv)){
    return 1;
  }else{
    return 0;
  }  
}

void CalcAdvPixImpact(AdvPix pix, Impact *imp, Impact *piximp){
  
  // Calculate impact coordinates in respect to the 
  // pixel coordinate system
  double u, v;
  
  u = imp->position.x - pix.sx;
  v = imp->position.y - pix.sy;
  
  // Fill piximp fields
  piximp->time = imp->time;
  piximp->energy = imp->energy;
  piximp->ph_id = imp->ph_id;
  piximp->src_id = imp->src_id;
  piximp->position.x = u;
  piximp->position.y = v;
}

int AdvImpactList(AdvDet *det, Impact *imp, long **pixindex, Impact **piximp){
  
  // Duplicate the impact but transform the coordinates into
  // the detector coordinate system
  Impact detimp;
  
  detimp.time = imp->time;
  detimp.energy = imp->energy;
  detimp.ph_id = imp->ph_id;
  detimp.src_id = imp->src_id;
  detimp.position.x = imp->position.x - det->sx;
  detimp.position.y = imp->position.y - det->sy;
  
  int nimpacts=0;
  
  // loop over all pixels and check for hit
  int ii;
  
  for(ii=0; ii<det->npix; ii++){
    if(CheckAdvPixImpact(det->pix[ii], &detimp)!=0){
      nimpacts++;
      *pixindex=(long*)realloc(*pixindex, nimpacts*sizeof(long));
      *pixindex[nimpacts-1]=(long)ii;
      *piximp=(Impact*)realloc(*piximp, nimpacts*sizeof(Impact));
      CalcAdvPixImpact(det->pix[ii], &detimp, &((*piximp)[nimpacts-1]));
    }
  }
  return nimpacts;
}

void parseAdvDetXML(AdvDet* const det, 
	       const char* const filename,
	       int* const status){

  headas_chat(5, "read advanced detector setup from XML file '%s' ...\n", filename);
		 
  // Read the XML data from the file.
  // Open the specified file.
  printf("Open file %s\n", filename);
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

  // Before expanding loops in the XML file, add the included code to it.
  expandIncludesXML(xmlbuffer, filename, status);
  CHECK_STATUS_VOID(*status);

  // Before acutally parsing the XML code, expand the loops and 
  // arithmetic operations in the XML description.
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
  struct XMLParseData xmlparsedata={
    .det  =det,
    .status=EXIT_SUCCESS
  };
  XML_SetUserData(parser, &xmlparsedata);

  // Set the handler functions.
  XML_SetElementHandler(parser, AdvDetXMLElementStart, AdvDetXMLElementEnd);

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
}

static void AdvDetXMLElementStart(void* parsedata, 
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
  
  // Check for advanced detector definition.
  if (!strcmp(Uelement, "PIXDETECTOR")) {
    // Determine npix
    char npix[MAXMSG];
    getXMLAttributeString(attr, "NPIX", npix);
    xmlparsedata->det->npix=atoi(npix);
    if(xmlparsedata->det->npix<1){
      SIXT_ERROR("Number of pixels in advanced detector description less than 1.");
      return;
    }
    xmlparsedata->det->pix=(AdvPix*)malloc(xmlparsedata->det->npix*sizeof(AdvPix));
    if(xmlparsedata->det->pix==NULL){
      SIXT_ERROR("Unable to allocate memory for advanced detector pixel array.");
      return;
    }
    xmlparsedata->det->sx=getXMLAttributeDouble(attr, "XOFF");
    xmlparsedata->det->sy=getXMLAttributeDouble(attr, "YOFF");    
  } else if (!strcmp(Uelement, "PIXEL")) {
    int posx, posy;
    posx=getXMLAttributeInt(attr, "POSX");
    posy=getXMLAttributeInt(attr, "POSY");
    xmlparsedata->det->pix[xmlparsedata->det->cpix].sx=(double)posx*getXMLAttributeDouble(attr, "DELX");
    xmlparsedata->det->pix[xmlparsedata->det->cpix].sy=(double)posy*getXMLAttributeDouble(attr, "DELY");
    xmlparsedata->det->pix[xmlparsedata->det->cpix].width=getXMLAttributeDouble(attr, "WIDTH");
    xmlparsedata->det->pix[xmlparsedata->det->cpix].height=getXMLAttributeDouble(attr, "HEIGHT");
    xmlparsedata->det->cpix++;
  } else {
    // Unknown tag, display warning.
    char msg[MAXMSG];
    sprintf(msg, "unknown XML tag: <%s>", el);
    SIXT_WARNING(msg);
  }
}

static void AdvDetXMLElementEnd(void* parsedata, const char* el) 
{
  struct XMLParseData* xmlparsedata=(struct XMLParseData*)parsedata;

  (void)el;

  // Check if an error has occurred previously.
  CHECK_STATUS_VOID(xmlparsedata->status);
}

AdvDet* loadAdvDet(const char* const filename,
		     int* const status)
{
  // Get a new and empty data structure.
  AdvDet* det=newAdvDet(status);
  CHECK_STATUS_RET(*status, det);

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
  CHECK_STATUS_RET(*status, det);

  // Split rootname into the file path and the file name.
  char* lastslash=strrchr(rootname, '/');
  if (NULL==lastslash) {
    det->filepath=(char*)malloc(sizeof(char));
    CHECK_NULL_RET(det->filepath, *status, 
		   "memory allocation for filepath failed", det);
    det->filename=(char*)malloc((strlen(rootname)+1)*sizeof(char));
    CHECK_NULL_RET(det->filename, *status, 
		   "memory allocation for filename failed", det);
    strcpy(det->filepath, "");
    strcpy(det->filename, rootname);
  } else {
    lastslash++;
    det->filename=(char*)malloc((strlen(lastslash)+1)*sizeof(char));
    CHECK_NULL_RET(det->filename, *status, 
		   "memory allocation for filename failed", det);
    strcpy(det->filename, lastslash);
      
    *lastslash='\0';
    det->filepath=(char*)malloc((strlen(rootname)+1)*sizeof(char));
    CHECK_NULL_RET(det->filepath, *status, 
		   "memory allocation for filepath failed", det);
    strcpy(det->filepath, rootname);
  }
  // END of storing the filename and filepath.


  // Read in the XML definition of the detector.
  parseAdvDetXML(det, filename, status);
  CHECK_STATUS_RET(*status, det);
  
  return(det);
}

