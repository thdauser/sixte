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
  module->xdim = 0.;
  for (ii=0; ii<module->nx; ii++) {
    // Determine the element with the maximum extension in x-direction 
    // (within the current column).
    float xmax=0.;
    long jj;
    for (jj=0; jj<module->ny; jj++) {
      if (module->element[ii+jj*module->nx]->xdim > xmax) {
	xmax = module->element[ii+jj*module->nx]->xdim;
      }
    }
    // Add the extension of the biggest element in this column.
    module->xdim += xmax;
  }

  // Loop over all rows.
  module->ydim = 0.;
  for (ii=0; ii<module->ny; ii++) {
    // Determine the element with the maximum extension in y-direction 
    // (within the current row).
    float ymax=0.;
    long jj;
    for (jj=0; jj<module->nx; jj++) {
      if (module->element[jj+ii*module->nx]->ydim > ymax) {
	ymax = module->element[jj+ii*module->nx]->ydim;
      }
    }
    // Add the extension of the biggest element in this row.
    module->ydim += ymax;
  }  
}


static void calcPanelXYDim(LADPanel* const panel)
{
  // XDIM and YDIM have to be calculated.
  long ii;

  // Loop over all columns.
  panel->xdim = 0.;
  for (ii=0; ii<panel->nx; ii++) {
    // Determine the module with the maximum extension in x-direction 
    // (within the current column).
    float xmax=0.;
    long jj;
    for (jj=0; jj<panel->ny; jj++) {
      // Calculate XDIM and YDIM for each child module.
      calcModuleXYDim(panel->module[ii+jj*panel->nx]);

      if (panel->module[ii+jj*panel->nx]->xdim > xmax) {
	xmax = panel->module[ii+jj*panel->nx]->xdim;
      }
    }
    // Add the extension of the biggest module in this column.
    panel->xdim += xmax;
  }

  // Loop over all rows.
  panel->ydim = 0.;
  for (ii=0; ii<panel->ny; ii++) {
    // Determine the module with the maximum extension in y-direction 
    // (within the current row).
    float ymax=0.;
    long jj;
    for (jj=0; jj<panel->nx; jj++) {
      if (panel->module[jj+ii*panel->nx]->ydim > ymax) {
	ymax = panel->module[jj+ii*panel->nx]->ydim;
      }
    }
    // Add the extension of the biggest module in this row.
    panel->ydim += ymax;
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
	
	headas_chat(5, "   element %ld has %ld anodes \n", 
		    lad->panel[ii]->module[jj]->element[kk]->id,
		    lad->panel[ii]->module[jj]->element[kk]->nanodes);
      }
      // END of loop over all elements.
    }
    // END of loop over all modules.


    // Calculate the dimensions of the panels from bottom up 
    // (elements -> modules -> panels).
    calcPanelXYDim(lad->panel[ii]);


    // Determine the geometric area of all panels and construct an ARF
    // from this value.
    float area=0.;
    for (ii=0; ii<lad->npanels; ii++) {
      area += lad->panel[ii]->xdim * lad->panel[ii]->ydim;
    }

    // Initialize an empty ARF data structure.
    lad->arf = getARF(status);
    CHECK_STATUS_VOID(*status);

    // Fill the ARF with data.
    lad->arf->NumberEnergyBins = 1;
    lad->arf->LowEnergy  = (float*)malloc(lad->arf->NumberEnergyBins*sizeof(float));
    CHECK_NULL_VOID(lad->arf->LowEnergy, *status, "memory allocation for ARF failed");
    lad->arf->HighEnergy = (float*)malloc(lad->arf->NumberEnergyBins*sizeof(float));
    CHECK_NULL_VOID(lad->arf->HighEnergy, *status, "memory allocation for ARF failed");
    lad->arf->EffArea    = (float*)malloc(lad->arf->NumberEnergyBins*sizeof(float));
    CHECK_NULL_VOID(lad->arf->EffArea, *status, "memory allocation for ARF failed");
    lad->arf->LowEnergy[0]     = 0.;
    lad->arf->HighEnergy[0]    = 1000.;
    lad->arf->EffArea[0]       = area;

  }
  // END of loop over all panels.
}


static void addString2XMLBuffer(struct XMLBuffer* const buffer, 
				const char* const string,
				int* const status)
{
  // Check if a valid buffer is specified.
  if (NULL==buffer) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("NULL pointer to XMLBuffer");
    return;
  }
    
  // Check if the buffer is empty.
  if (NULL==buffer->text) {
    // Allocate memory for the first chunk of bytes.
    buffer->text=(char*)malloc((MAXMSG+1)*sizeof(char));
    if (NULL==buffer->text) {
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for XMLBuffer failed");
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
      SIXT_ERROR("memory allocation for XMLBuffer failed");
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
    SIXT_ERROR("memory allocation for XMLBuffer failed");
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
    SIXT_ERROR("memory allocation for XMLBuffer failed");
    return(buffer);
  }

  buffer->text=NULL;
  buffer->maxlength=0;

  return(buffer);
}


static void freeXMLBuffer(struct XMLBuffer** const buffer)
{
  if (NULL!=*buffer) {
    if (NULL!=(*buffer)->text) {
      free((*buffer)->text);
    }
    free(*buffer);
    *buffer=NULL;
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
  lad->panel[lad->npanels-1] = panel;
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
  panel->module = 
    (LADModule**)realloc(panel->module, (panel->nmodules+1)*sizeof(LADModule*));
  CHECK_NULL_VOID(panel->module, *status, "memory allocation for new LADModule failed");
  panel->nmodules++;

  // Append the new module to the LAD.
  panel->module[panel->nmodules-1] = module;
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
  module->element = 
    (LADElement**)realloc(module->element, (module->nelements+1)*sizeof(LADElement*));
  CHECK_NULL_VOID(module->element, *status, "memory allocation for new LADElement failed");
  module->nelements++;

  // Append the new element to the LAD.
  module->element[module->nelements-1] = element;
}


static void XMLElementStart(void* parsedata, const char* el, const char** attr) 
{
  struct XMLParseData* xmlparsedata = (struct XMLParseData*)parsedata;
  char Uelement[MAXMSG];   // Upper case version of XML element

  // Check if an error has occurred previously.
  CHECK_STATUS_VOID(xmlparsedata->status);

  // Convert the element to an upper case string.
  strcpy(Uelement, el);
  strtoupper(Uelement);

  // Check for the different kinds of allowed elements.
  if (!strcmp(Uelement, "PANEL")) {

    // Determine the ID of the new panel.
    char buffer[MAXMSG]; // String buffer.
    getAttribute(attr, "ID", buffer);
    long id = atol(buffer);
    getAttribute(attr, "NX", buffer);
    long nx = atol(buffer);
    getAttribute(attr, "NY", buffer);
    long ny = atol(buffer);

    // Create a new Panel.
    LADPanel* panel = newLADPanel(&xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);
    
    // Set the properties.
    panel->id = id;
    panel->nx = nx;
    panel->ny = ny;
    
    // Add the new panel to the LAD.
    addPanel2LAD(xmlparsedata->lad, panel, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);
    
    // Store the pointer to the currently open panel.
    xmlparsedata->panel = panel;

  } else if (!strcmp(Uelement, "MODULE")) {

    // Determine the ID of the new module.
    char buffer[MAXMSG]; // String buffer.
    getAttribute(attr, "ID", buffer);
    long id = atol(buffer);
    getAttribute(attr, "NX", buffer);
    long nx = atol(buffer);
    getAttribute(attr, "NY", buffer);
    long ny = atol(buffer);

    // Create a new Module.
    LADModule* module = newLADModule(&xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);
    
    // Set the properties.
    module->id = id;
    module->nx = nx;
    module->ny = ny;

    // Add the new module to the LAD.
    addModule2Panel(xmlparsedata->panel, module, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);
    
    // Store the pointer to the currently open module.
    xmlparsedata->module = module;

  } else if (!strcmp(Uelement, "ELEMENT")) {

    // Determine the ID of the new element.
    char buffer[MAXMSG]; // String buffer.
    getAttribute(attr, "ID", buffer);
    long id = atol(buffer);
    getAttribute(attr, "XDIM", buffer);
    float xdim = (float)atof(buffer);
    getAttribute(attr, "YDIM", buffer);
    float ydim = (float)atof(buffer);
    getAttribute(attr, "NANODES", buffer);
    long nanodes = atol(buffer);
    getAttribute(attr, "ANODEPITCH", buffer);
    float anodepitch = (float)atof(buffer);
    
    // Create a new Element.
    LADElement* element = newLADElement(&xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);
    
    // Set the properties.
    element->id = id;
    element->xdim = xdim;
    element->ydim = ydim;
    element->nanodes = nanodes;
    element->anodepitch = anodepitch;
    
    // Add the new element to the LAD.
    addElement2Module(xmlparsedata->module, element, &xmlparsedata->status);
    CHECK_STATUS_VOID(xmlparsedata->status);
    
    // Store the pointer to the currently open element.
    xmlparsedata->element = element;

  } else {
    xmlparsedata->status = EXIT_FAILURE;
    SIXT_ERROR("unknown XML element");
  }
}


static void XMLElementEnd(void* parsedata, const char* el) 
{
  struct XMLParseData* xmlparsedata = (struct XMLParseData*)parsedata;
  char Uelement[MAXMSG]; // Upper case version of XML element

  // Check if an error has occurred previously.
  CHECK_STATUS_VOID(xmlparsedata->status);

  // Convert the element to an upper case string.
  strcpy(Uelement, el);
  strtoupper(Uelement);

  // Check, whether a panel, module, or element has ben terminated.
  if (!strcmp(Uelement, "PANEL")) {
    xmlparsedata->panel = NULL;
  } else if (!strcmp(Uelement, "MODULE")) {
    xmlparsedata->module = NULL;
  } else if (!strcmp(Uelement, "ELEMENT")) {
    xmlparsedata->element = NULL;
  }

  return;
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
      SIXT_ERROR("memory allocation for string buffer in XML "
		 "pre-parsing failed");
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
    SIXT_ERROR("XML element string too long");
    return;
  }
  addString2XMLBuffer(output, buffer, &mydata->status);
  CHECK_STATUS_VOID(mydata->status);

  int ii=0;
  while(attr[ii]) {
    if (sprintf(buffer, " %s=\"%s\"", attr[ii], attr[ii+1]) >= MAXMSG) {
      mydata->status=EXIT_FAILURE;
      SIXT_ERROR("XML element string too long");
      return;
    }
    addString2XMLBuffer(output, buffer, &mydata->status);
    CHECK_STATUS_VOID(mydata->status);
    
    ii+=2;
  }

  addString2XMLBuffer(output, ">", &mydata->status);
  CHECK_STATUS_VOID(mydata->status);
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
      SIXT_ERROR("memory allocation for string buffer in XML "
		     "pre-parsing failed");
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
	CHECK_STATUS_VOID(mydata->status);
	copyXMLBuffer(replacedBuffer, mydata->loop_buffer, &mydata->status);
	CHECK_STATUS_VOID(mydata->status);

	// Replace $variables by integer values.
	if (strlen(replacedBuffer->text)>0) {
	  char stringvalue[MAXMSG];
	  sprintf(stringvalue, "%d", ii);
	  replaceInXMLBuffer(replacedBuffer, mydata->loop_variable,
			     stringvalue, &mydata->status);
	  CHECK_STATUS_VOID(mydata->status);
	}

	// Add the loop content to the output buffer.
	addString2XMLBuffer(mydata->output_buffer, replacedBuffer->text, 
			    &mydata->status);
	CHECK_STATUS_VOID(mydata->status);

	freeXMLBuffer(&replacedBuffer);
      }
      // Clear the loop buffer.
      freeXMLBuffer(&mydata->loop_buffer);
      mydata->loop_buffer=newXMLBuffer(&mydata->status);
      CHECK_STATUS_VOID(mydata->status);
      
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
    SIXT_ERROR("XML string element too long");
    return;
  }
  addString2XMLBuffer(output, buffer, &mydata->status);
  CHECK_STATUS_VOID(mydata->status);
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
      SIXT_ERROR("could not allocate memory for XML parser");
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
      sprintf(msg, "parsing XML code failed: %s", 
	      XML_ErrorString(XML_GetErrorCode(parser)));
      SIXT_ERROR(msg);
      return;
    }
    // Check for errors.
    if (EXIT_SUCCESS!=data.status) {
      *status = data.status;
      return;
    }

    // Copy the output XMLBuffer to the input XMLBuffer ...
    copyXMLBuffer(buffer, data.output_buffer, status);
    CHECK_STATUS_VOID(*status);
    // ... and release allocated memory.
    freeXMLBuffer(&data.output_buffer);
    freeXMLBuffer(&data.loop_buffer);

    XML_ParserFree(parser);

  } while (data.further_loops>0);

  // Replace arithmetic +/- expressions.
  execArithmeticOpsInXMLBuffer(buffer, status);
}


LAD* getLADfromXML(const char* const filename, 
		   int* const status)
{
  LAD* lad=NULL;

  // Allocate memory for a new LAD data structure.
  lad = newLAD(status);
  CHECK_STATUS_RET(*status, lad);


  // Read the data from the XML file.
  // Open the XML file.
  FILE* xmlfile = fopen(filename, "r");
  if (NULL==xmlfile) {
    *status = EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "failed opening LAD XML definition "
	    "file '%s' for read access!\n", filename);
    SIXT_ERROR(msg);
    return(lad);
  }

  // The data is read from the XML file and stored in xmlbuffer
  // without any modifications.
  struct XMLBuffer* xmlbuffer = newXMLBuffer(status);
  CHECK_STATUS_RET(*status, lad);

  // Input buffer with an additional byte at the end for the 
  // termination of the string.
  char buffer[MAXMSG+1];
  // Number of chars in buffer.
  int len;

  // Read all data from the file.
  do {
    // Get a piece of input into the buffer.
    len = fread(buffer, 1, MAXMSG, xmlfile);
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
  XML_Parser parser = XML_ParserCreate(NULL);
  if (NULL==parser) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("could not allocate memory for XML parser");
    return(lad);
  }

  // Set data that is passed to the handler functions.
  struct XMLParseData xmlparsedata = {
    .lad    = lad,
    .panel  = NULL,
    .module = NULL,
    .element= NULL,
    .status = EXIT_SUCCESS
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
  CHECK_STATUS_RET(xmlparsedata.status, lad);
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

