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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "lad.h"


////////////////////////////////////////////////////////////////////
// Local data type declarations.
////////////////////////////////////////////////////////////////////


/** Data structure given to the XML handler to transfer data. */
struct XMLParseData {
  LAD* lad;
  LADPanel* panel;
  LADModule* module;
  LADElement* element;

  int status;
};


////////////////////////////////////////////////////////////////////
// Static function declarations.
////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////
// Program Code.
////////////////////////////////////////////////////////////////////


static void calcModuleXYDim(LADModule* const module)
{
  // XDIM and YDIM have to be calculated.
  long ii;

  // Loop over all columns.
  module->xdim=0.;
  for (ii=0; ii<module->nx; ii++) {
    // Determine the element with the maximum extension in x-direction
    // (within the current column).
    float xmax=0.;
    long jj;
    for (jj=0; jj<module->ny; jj++) {
      if (module->element[ii+jj*module->nx]->xdim > xmax) {
	xmax=module->element[ii+jj*module->nx]->xdim;
      }
    }
    // Add the extension of the biggest element in this column.
    module->xdim+=xmax;
  }

  // Loop over all rows.
  module->ydim=0.;
  for (ii=0; ii<module->ny; ii++) {
    // Determine the element with the maximum extension in y-direction
    // (within the current row).
    float ymax=0.;
    long jj;
    for (jj=0; jj<module->nx; jj++) {
      if (module->element[jj+ii*module->nx]->ydim > ymax) {
	ymax=module->element[jj+ii*module->nx]->ydim;
      }
    }
    // Add the extension of the biggest element in this row.
    module->ydim+=ymax;
  }
}


static void calcPanelXYDim(LADPanel* const panel)
{
  // XDIM and YDIM have to be calculated.
  long ii;

  // Loop over all columns.
  panel->xdim=0.;
  for (ii=0; ii<panel->nx; ii++) {
    // Determine the module with the maximum extension in x-direction
    // (within the current column).
    float xmax=0.;
    long jj;
    for (jj=0; jj<panel->ny; jj++) {
      // Calculate XDIM and YDIM for each child module.
      calcModuleXYDim(panel->module[ii+jj*panel->nx]);

      if (panel->module[ii+jj*panel->nx]->xdim > xmax) {
	xmax=panel->module[ii+jj*panel->nx]->xdim;
      }
    }
    // Add the extension of the biggest module in this column.
    panel->xdim+=xmax;
  }

  // Loop over all rows.
  panel->ydim=0.;
  for (ii=0; ii<panel->ny; ii++) {
    // Determine the module with the maximum extension in y-direction
    // (within the current row).
    float ymax=0.;
    long jj;
    for (jj=0; jj<panel->nx; jj++) {
      if (panel->module[jj+ii*panel->nx]->ydim > ymax) {
	ymax=panel->module[jj+ii*panel->nx]->ydim;
      }
    }
    // Add the extension of the biggest module in this row.
    panel->ydim+=ymax;
  }
}


static void checkLADConsistency(LAD* const lad, int* const status)
{
  // Check if the the LAD exists.
  CHECK_NULL_VOID(lad, *status, "NULL pointer to LAD data structure");

  headas_chat(5, "LAD detector\n");

  // Check if the LAD contains any panels.
  CHECK_NULL_VOID(lad->panel, *status, "LAD contains no panels");
  // Loop over all panels.
  long ii;
  for (ii=0; ii<lad->npanels; ii++) {
    // Check if the the panel exists.
    CHECK_NULL_VOID(lad->panel[ii], *status, "panel is not defined");

    // Check if the panel contains any modules.
    CHECK_NULL_VOID(lad->panel[ii]->module, *status,
		    "panel contains no modules");

    // Check if the number of modules is consistent with the
    // product of nx times ny.
    if (lad->panel[ii]->nmodules!=lad->panel[ii]->nx*lad->panel[ii]->ny) {
      *status = EXIT_FAILURE;
      SIXT_ERROR("inconsistent number of modules in a panel");
      return;
    }

    headas_chat(5, " panel %ld contains %ld modules\n",
		lad->panel[ii]->id, lad->panel[ii]->nmodules);

    // Loop over all modules.
    long jj;
    for (jj=0; jj<lad->panel[ii]->nmodules; jj++) {
      // Check if the the module exists.
      CHECK_NULL_VOID(lad->panel[ii]->module[jj], *status,
		      "module is not defined");

      // Check if the module contains any elements.
      CHECK_NULL_VOID(lad->panel[ii]->module[jj]->element, *status,
		      "module contains no elements");

      // Check if the number of elements is consistent with the
      // product of nx times ny.
      if (lad->panel[ii]->module[jj]->nelements!=
	  lad->panel[ii]->module[jj]->nx*lad->panel[ii]->module[jj]->ny) {
	*status = EXIT_FAILURE;
	SIXT_ERROR("inconsistent number of elements in a module");
	return;
      }

      headas_chat(5, "  module %ld contains %ld elements\n",
		  lad->panel[ii]->module[jj]->id,
		  lad->panel[ii]->module[jj]->nelements);

      // Loop over all elements.
      long kk;
      for (kk=0; kk<lad->panel[ii]->module[jj]->nelements; kk++) {
	// Check if the the element exists.
	CHECK_NULL_VOID(lad->panel[ii]->module[jj]->element[kk],
			*status, "element is not defined");

	// Check if the element contains any anodes.
	if (0==lad->panel[ii]->module[jj]->element[kk]->nanodes) {
	  SIXT_ERROR("element contains no anodes");
	  *status=EXIT_FAILURE;
	  return;
	}

	// Check if the number of anodes is even.
	if (lad->panel[ii]->module[jj]->element[kk]->nanodes % 2 == 1) {
	  SIXT_ERROR("number of anodes must be even");
	  *status=EXIT_FAILURE;
	  return;
	}

	// Set up the FEE ASICs.
	// Determine the required number of ASICs.
	if (lad->panel[ii]->module[jj]->element[kk]->nanodes/2 %
	    lad->asic_channels != 0) {
	  SIXT_ERROR("number of anodes must be a multiple of ASIC channels");
	  *status=EXIT_FAILURE;
	  return;
	}
	lad->panel[ii]->module[jj]->element[kk]->nasics=
	  (int)(lad->panel[ii]->module[jj]->element[kk]->nanodes/
		lad->asic_channels);
	// Allocate memory for the read-out times of the ASICs.
	lad->panel[ii]->module[jj]->element[kk]->asic_readout_time=
	  (double*)malloc(lad->panel[ii]->module[jj]->element[kk]->nasics*
			  sizeof(double));
	CHECK_NULL_VOID(lad->panel[ii]->module[jj]->element[kk]->
			asic_readout_time, *status,
			"cannot allocate memory for ASIC read-out times");
	int ll;
	for (ll=0; ll<lad->panel[ii]->module[jj]->element[kk]->nasics; ll++) {
	  lad->panel[ii]->module[jj]->element[kk]->asic_readout_time[ll]=0.;
	}

	// Allocate memory for the dead times of the individual ASICs.
	lad->panel[ii]->module[jj]->element[kk]->asic_deadtime=
	  (double*)malloc(lad->panel[ii]->module[jj]->element[kk]->nasics*
			  sizeof(double));
	CHECK_NULL_VOID(lad->panel[ii]->module[jj]->element[kk]->
			asic_deadtime, *status,
			"cannot allocate memory for ASIC dead times");
	for (ll=0; ll<lad->panel[ii]->module[jj]->element[kk]->nasics; ll++) {
	  double grand1, grand2;
	  sixt_get_gauss_random_numbers(&grand1, &grand2, status);
	  CHECK_STATUS_VOID(*status);
	  lad->panel[ii]->module[jj]->element[kk]->asic_deadtime[ll]=
	    lad->deadtime + lad->edeadtime*grand1;
	}

	headas_chat(5, "   element %ld has %ld anodes and %d ASICs \n",
		    lad->panel[ii]->module[jj]->element[kk]->id,
		    lad->panel[ii]->module[jj]->element[kk]->nanodes,
		    lad->panel[ii]->module[jj]->element[kk]->nasics);
      }
      // END of loop over all elements.
    }
    // END of loop over all modules.


    // Calculate the dimensions of the panels from bottom up
    // (elements -> modules -> panels).
    calcPanelXYDim(lad->panel[ii]);

  }
  // END of loop over all panels.
}


static void addPanel2LAD(LAD* const lad,
			 LADPanel* const panel,
			 int* const status)
{
  // Check if the LAD is defined.
  CHECK_NULL_VOID(lad, *status, "NULL pointer to LAD data structure");

  // Check if the panel is defined.
  CHECK_NULL_VOID(panel, *status, "NULL pointer to LADPanel data structure");

  // Extend the LAD panel array.
  lad->panel =
    (LADPanel**)realloc(lad->panel, (lad->npanels+1)*sizeof(LADPanel*));
  CHECK_NULL_VOID(lad->panel, *status, "memory allocation for new LADPanel failed");
  lad->npanels++;

  // Append the new panel to the LAD.
  lad->panel[lad->npanels-1]=panel;
}


static void addModule2Panel(LADPanel* const panel,
			    LADModule* const module,
			    int* const status)
{
  // Check if the panel is defined.
  CHECK_NULL_VOID(panel, *status, "NULL pointer to LADPanel data structure");

  // Check if the module is defined.
  CHECK_NULL_VOID(module, *status, "NULL pointer to LADModule data structure");

  // Extend the LAD module array.
  panel->module=
    (LADModule**)realloc(panel->module, (panel->nmodules+1)*sizeof(LADModule*));
  CHECK_NULL_VOID(panel->module, *status,
		  "memory allocation for new LADModule failed");
  panel->nmodules++;

  // Append the new module to the LAD.
  panel->module[panel->nmodules-1]=module;
}


static void addElement2Module(LADModule* const module,
			      LADElement* const element,
			      int* const status)
{
  // Check if the module is defined.
  CHECK_NULL_VOID(module, *status, "NULL pointer to LADModule data structure");

  // Check if the element is defined.
  CHECK_NULL_VOID(element, *status, "NULL pointer to LADElement data structure");

  // Extend the LAD element array.
  module->element=
    (LADElement**)realloc(module->element,
			  (module->nelements+1)*sizeof(LADElement*));
  CHECK_NULL_VOID(module->element, *status,
		  "memory allocation for new LADElement failed");
  module->nelements++;

  // Append the new element to the LAD.
  module->element[module->nelements-1]=element;
}


static void XMLElementStart(void* parsedata, const char* el, const char** attr)
{
  struct XMLParseData* xmlparsedata=(struct XMLParseData*)parsedata;
  char Uelement[MAXMSG]; // Upper case version of XML element.

  // Check if an error has occurred previously.
  CHECK_STATUS_VOID(xmlparsedata->status);

  // Convert the element to an upper case string.
  strcpy(Uelement, el);
  strtoupper(Uelement);

  // Check for the different kinds of allowed elements.
  if (!strcmp(Uelement, "PANEL")) {

    // Create a new Panel.
    LADPanel* panel=newLADPanel(&xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

    // Set the properties.
    panel->id=getXMLAttributeLong(attr, "ID");
    panel->nx=getXMLAttributeLong(attr, "NX");
    panel->ny=getXMLAttributeLong(attr, "NY");

    // Add the new panel to the LAD.
    addPanel2LAD(xmlparsedata->lad, panel, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

    // Store the pointer to the currently open panel.
    xmlparsedata->panel=panel;

  } else if (!strcmp(Uelement, "MODULE")) {

    // Create a new Module.
    LADModule* module=newLADModule(&xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

    // Set the properties.
    module->id=getXMLAttributeLong(attr, "ID");
    module->nx=getXMLAttributeLong(attr, "NX");
    module->ny=getXMLAttributeLong(attr, "NY");

    // Add the new module to the LAD.
    addModule2Panel(xmlparsedata->panel, module, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

    // Store the pointer to the currently open module.
    xmlparsedata->module=module;

  } else if (!strcmp(Uelement, "ELEMENT")) {

    // Create a new Element.
    LADElement* element=newLADElement(&xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

    // Set the properties.
    element->id=getXMLAttributeLong(attr, "ID");
    element->xdim=getXMLAttributeFloat(attr, "XDIM");
    element->ydim=getXMLAttributeFloat(attr, "YDIM");
    element->xborder=getXMLAttributeFloat(attr, "XBORDER");
    element->yborder=getXMLAttributeFloat(attr, "YBORDER");
    element->nanodes=getXMLAttributeLong(attr, "NANODES");

    // Add the new element to the LAD.
    addElement2Module(xmlparsedata->module, element, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

    // Store the pointer to the currently open element.
    xmlparsedata->element=element;

  } else if (!strcmp(Uelement, "FOV")) {

    // Determine the diameter of the FOV.
    xmlparsedata->lad->fov_diameter=
      getXMLAttributeFloat(attr, "DIAMETER")*M_PI/180.;

  } else if (!strcmp(Uelement, "TEMPERATURE")) {

    // Determine the value of the temperature.
    xmlparsedata->lad->temperature=getXMLAttributeFloat(attr, "VALUE");

  } else if (!strcmp(Uelement, "EFIELD")) {

    // Determine the electric field.
    xmlparsedata->lad->efield=getXMLAttributeFloat(attr, "VALUE");

  } else if (!strcmp(Uelement, "MOBILITY")) {

    // Determine the mobility.
    xmlparsedata->lad->mobility=getXMLAttributeFloat(attr, "VALUE");

  } else if (!strcmp(Uelement, "ASIC")) {

    // Determine the dead time.
    xmlparsedata->lad->deadtime=
      getXMLAttributeFloat(attr, "DEADTIME");
    // Determine the error on the dead time
    xmlparsedata->lad->edeadtime=
      getXMLAttributeFloat(attr, "EDEADTIME");
    // Determine the length of the coincidence time window.
    xmlparsedata->lad->coincidencetime=
      getXMLAttributeFloat(attr, "COINCIDENCETIME");
    // Determine the number of input channels per ASIC.
    xmlparsedata->lad->asic_channels=
      getXMLAttributeInt(attr, "CHANNELS");

  } else if (!strcmp(Uelement, "THRESHOLD")) {

    // Determine the lower and the upper read-out threshold.
    float threshold_readout_lo_keV=
      getXMLAttributeFloat(attr, "READOUT_LO_KEV");
    float threshold_readout_up_keV=
      getXMLAttributeFloat(attr, "READOUT_UP_KEV");

    if (0.<threshold_readout_lo_keV) {
      xmlparsedata->lad->threshold_readout_lo_keV=
	(float*)malloc(sizeof(float));
      CHECK_NULL_VOID(xmlparsedata->lad->threshold_readout_lo_keV,
		      xmlparsedata->status,
		      "memory allocation for threshold failed");
      *(xmlparsedata->lad->threshold_readout_lo_keV)=threshold_readout_lo_keV;
    }

    if (0.<threshold_readout_up_keV) {
      xmlparsedata->lad->threshold_readout_up_keV=
	(float*)malloc(sizeof(float));
      CHECK_NULL_VOID(xmlparsedata->lad->threshold_readout_up_keV,
		      xmlparsedata->status,
		      "memory allocation for threshold failed");
      *(xmlparsedata->lad->threshold_readout_up_keV)=threshold_readout_up_keV;
    }

  } else if (!strcmp(Uelement, "ARF")) {

    // Check if the ARF has been defined previously.
    if (NULL!=xmlparsedata->lad->arf) {
      xmlparsedata->status=EXIT_FAILURE;
      SIXT_ERROR("ARF already defined (cannot be loaded twice)");
      return;
    }

    // Determine the file name of the instrument ARF.
    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);

    // Check if a file name has been specified.
    if (strlen(filename)==0) {
      xmlparsedata->status=EXIT_FAILURE;
      SIXT_ERROR("no file specified for ARF");
      return;
    }

    // Store the file name of the ARF.
    xmlparsedata->lad->arf_filename=
      (char*)malloc((strlen(filename)+1)*sizeof(char));
    CHECK_NULL_VOID(xmlparsedata->lad->arf_filename,
		    xmlparsedata->status,
		    "memory allocation for ARF file name failed");
    strcpy(xmlparsedata->lad->arf_filename, filename);

    // Load the ARF.
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->lad->filepath);
    strcat(filepathname, filename);
    xmlparsedata->lad->arf=loadARF(filepathname, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

  } else if (!strcmp(Uelement, "RMF")) {

    // Check if the RMF has been defined previously.
    if (NULL!=xmlparsedata->lad->rmf) {
      SIXT_ERROR("RMF already defined (cannot be loaded twice)");
      xmlparsedata->status=EXIT_FAILURE;
      return;
    }

    // Determine the file name of the RMF.
    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);

    // Check if a file name has been specified.
    if (strlen(filename)==0) {
      SIXT_ERROR("no file specified for RMF");
      xmlparsedata->status=EXIT_FAILURE;
      return;
    }

    // Store the file name of the RMF.
    xmlparsedata->lad->rmf_filename=
      (char*)malloc((strlen(filename)+1)*sizeof(char));
    CHECK_NULL_VOID(xmlparsedata->lad->rmf_filename,
		    xmlparsedata->status,
		    "memory allocation for RMF file name failed");
    strcpy(xmlparsedata->lad->rmf_filename, filename);

    // Load the RMF.
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->lad->filepath);
    strcat(filepathname, filename);
    xmlparsedata->lad->rmf=loadNormalizedRMF(filepathname, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

  } else if (!strcmp(Uelement, "RSP")) {

    // Check if the ARF or RMF have been defined previously.
    if ((NULL!=xmlparsedata->lad->arf)||(NULL!=xmlparsedata->lad->rmf)) {
      SIXT_ERROR("ARF or RMF already defined (cannot be loaded twice)");
      xmlparsedata->status=EXIT_FAILURE;
      return;
    }

    // Determine the file name of the RSP.
    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);

    // Check if a file name has been specified.
    if (strlen(filename)==0) {
      SIXT_ERROR("no file specified for RSP");
      xmlparsedata->status=EXIT_FAILURE;
      return;
    }

    // Store the file name of the RSP.
    xmlparsedata->lad->arf_filename=
      (char*)malloc((strlen(filename)+1)*sizeof(char));
    CHECK_NULL_VOID(xmlparsedata->lad->arf_filename,
		    xmlparsedata->status,
		    "memory allocation for ARF file name failed");
    strcpy(xmlparsedata->lad->arf_filename, filename);
    xmlparsedata->lad->rmf_filename=
      (char*)malloc((strlen(filename)+1)*sizeof(char));
    CHECK_NULL_VOID(xmlparsedata->lad->rmf_filename,
		    xmlparsedata->status,
		    "memory allocation for RMF file name failed");
    strcpy(xmlparsedata->lad->rmf_filename, filename);

    // Load the RSP.
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->lad->filepath);
    strcat(filepathname, filename);
    loadArfRmfFromRsp(filepathname,
		      &xmlparsedata->lad->arf,
		      &xmlparsedata->lad->rmf,
		      &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

  } else if (!strcmp(Uelement, "BACKGROUND")) {

    // Determine the file name of the background catalog.
    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);

    // Check if a file name has been specified.
    if (strlen(filename)==0) {
      SIXT_ERROR("no file specified for detector background");
      xmlparsedata->status=EXIT_FAILURE;
    }

    // Open the background SIMPUT catalog file.
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->lad->filepath);
    strcat(filepathname, filename);
    xmlparsedata->lad->bkgctlg=
      openSimputCtlg(filepathname, READONLY,
		     0, 0, 0, 0, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

  } else if (!strcmp(Uelement, "VIGNETTING")) {

    // Determine the file name of the collimator vignetting function.
    char filename[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", filename);

    // Check if a file name has been specified.
    if (strlen(filename)==0) {
      SIXT_ERROR("no file specified for vignetting");
      xmlparsedata->status=EXIT_FAILURE;
    }

    // Load the vignetting function.
    char filepathname[MAXFILENAME];
    strcpy(filepathname, xmlparsedata->lad->filepath);
    strcat(filepathname, filename);
    xmlparsedata->lad->vignetting=
      newVignetting(filepathname, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);

  } else {
    char msg[MAXMSG];
    sprintf(msg, "unknown XML element '%s'", Uelement);
    SIXT_WARNING(msg);
  }
}


static void XMLElementEnd(void* parsedata, const char* el)
{
  struct XMLParseData* xmlparsedata=(struct XMLParseData*)parsedata;
  char Uelement[MAXMSG]; // Upper case version of XML element

  // Check if an error has occurred previously.
  CHECK_STATUS_VOID(xmlparsedata->status);

  // Convert the element to an upper case string.
  strcpy(Uelement, el);
  strtoupper(Uelement);

  // Check, whether a panel, module, or element has ben terminated.
  if (!strcmp(Uelement, "PANEL")) {
    xmlparsedata->panel=NULL;
  } else if (!strcmp(Uelement, "MODULE")) {
    xmlparsedata->module=NULL;
  } else if (!strcmp(Uelement, "ELEMENT")) {
    xmlparsedata->element=NULL;
  }

  return;
}


LAD* getLADfromXML(const char* const filename,
		   int* const status)
{
  LAD* lad=NULL;

  // Allocate memory for a new LAD data structure.
  lad=newLAD(status);
  CHECK_STATUS_RET(*status, lad);


  // Split the reference to the XML LAD definition file
  // into file path and name. This has to be done before
  // calling the parser routine for the XML file.
  char filename2[MAXFILENAME];
  char rootname[MAXFILENAME];
  // Make a local copy of the filename variable in order to avoid
  // compiler warnings due to discarded const qualifier at the
  // subsequent function call.
  strcpy(filename2, filename);
  fits_parse_rootname(filename2, rootname, status);
  CHECK_STATUS_RET(*status, lad);

  // Split rootname into the file path and the file name.
  char* lastslash = strrchr(rootname, '/');
  if (NULL==lastslash) {
    lad->filepath=(char*)malloc(sizeof(char));
    CHECK_NULL_RET(lad->filepath, *status,
		   "memory allocation for filepath failed", lad);
    lad->filename=(char*)malloc((strlen(rootname)+1)*sizeof(char));
    CHECK_NULL_RET(lad->filename, *status,
		   "memory allocation for filename failed", lad);
    strcpy(lad->filepath, "");
    strcpy(lad->filename, rootname);
  } else {
    lastslash++;
    lad->filename=(char*)malloc((strlen(lastslash)+1)*sizeof(char));
    CHECK_NULL_RET(lad->filename, *status,
		   "memory allocation for filename failed", lad);
    strcpy(lad->filename, lastslash);

    *lastslash='\0';
    lad->filepath=(char*)malloc((strlen(rootname)+1)*sizeof(char));
    CHECK_NULL_RET(lad->filepath, *status,
		   "memory allocation for filepath failed", lad);
    strcpy(lad->filepath, rootname);
  }
  // END of storing the filename and filepath.


  // Read the data from the XML file.
  // Open the XML file.
  FILE* xmlfile=fopen(filename, "r");
  if (NULL==xmlfile) {
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "failed opening LAD XML definition "
	    "file '%s' for read access!\n", filename);
    SIXT_ERROR(msg);
    return(lad);
  }

  // The data is read from the XML file and stored in xmlbuffer
  // without any modifications.
  struct XMLBuffer* xmlbuffer=newXMLBuffer(status);
  CHECK_STATUS_RET(*status, lad);

  // Input buffer with an additional byte at the end for the
  // termination of the string.
  char buffer[MAXMSG+1];
  // Number of chars in buffer.
  int len;

  // Read all data from the file.
  do {
    // Get a piece of input into the buffer.
    len=fread(buffer, 1, MAXMSG, xmlfile);
    buffer[len]='\0'; // Terminate the string.
    addString2XMLBuffer(xmlbuffer, buffer, status);
    CHECK_STATUS_RET(*status, lad);
  } while (!feof(xmlfile));

  // Close the file handler to the XML file.
  fclose(xmlfile);
  // END of reading the detector definition from the XML file.


  // Consistency check of XML code ???
  // TODO


  // Preprocess the XML code (expand loops, perform mathematical
  // operations).
  // Before acutally parsing the XML code, expand the loops and
  // arithmetic operations in the XML description.
  // The expansion algorithm repeatetly scans the XML code and
  // searches for loop tags. It replaces the loop tags by repeating
  // the contained XML code.
  expandXML(xmlbuffer, status);
  CHECK_STATUS_RET(*status, lad);


  // Iteratively parse the XML code and construct the LAD data
  // structure.
  // Parse XML code in the xmlbuffer using the expat library.
  // Get an XML_Parser object.
  XML_Parser parser=XML_ParserCreate(NULL);
  if (NULL==parser) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("could not allocate memory for XML parser");
    return(lad);
  }

  // Set data that is passed to the handler functions.
  struct XMLParseData xmlparsedata = {
    .lad    =lad,
    .panel  =NULL,
    .module =NULL,
    .element=NULL,
    .status =EXIT_SUCCESS
  };
  XML_SetUserData(parser, &xmlparsedata);

  // Set the handler functions.
  XML_SetElementHandler(parser, XMLElementStart, XMLElementEnd);

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
    return(lad);
  }
  // Check for errors.
  *status=xmlparsedata.status;
  CHECK_STATUS_RET(*status, lad);

  // Release memory.
  XML_ParserFree(parser);

  // Remove the XML string buffer.
  freeXMLBuffer(&xmlbuffer);

  // Iteratively go through the LAD data structure and set properties
  // of parent and child elements. Perform consistency check of the
  // LAD detector.
  checkLADConsistency(lad, status);
  CHECK_STATUS_RET(*status, lad);

  return(lad);
}
