#include "xmlbuffer.h"


void addString2XMLBuffer(struct XMLBuffer* const buffer, 
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
    int new_length=strlen(buffer->text) + strlen(string);
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
  destination->text=(char*)realloc(destination->text,
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


struct XMLBuffer* newXMLBuffer(int* const status)
{
  struct XMLBuffer* buffer=
    (struct XMLBuffer*)malloc(sizeof(struct XMLBuffer));
  if (NULL==buffer) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for XMLBuffer failed");
    return(buffer);
  }

  buffer->text=NULL;
  buffer->maxlength=0;

  return(buffer);
}


void freeXMLBuffer(struct XMLBuffer** const buffer)
{
  if (NULL!=*buffer) {
    if (NULL!=(*buffer)->text) {
      free((*buffer)->text);
    }
    free(*buffer);
    *buffer=NULL;
  }
}


static void expandXMLElementStart(void* data, const char* el, 
				  const char** attr)
{
  struct XMLPreParseData* mydata=(struct XMLPreParseData*)data;

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
	  mydata->loop_start=atoi(attr[ii+1]);
	} else if (!strcmp(Uattribute, "END")) {
	  mydata->loop_end=atoi(attr[ii+1]);
	} else if (!strcmp(Uattribute, "INCREMENT")) {
	  mydata->loop_increment=atoi(attr[ii+1]);
	} else if (!strcmp(Uattribute, "VARIABLE")) {
	  strcpy(mydata->loop_variable, attr[ii+1]);
	}

	ii+=2;
      }
      // END of loop over all attributes.
      
      // Check if parameters are set to valid values.
      if (((mydata->loop_end-mydata->loop_start)*mydata->loop_increment<0) ||
	  (0==mydata->loop_increment)){
	mydata->status=EXIT_FAILURE;
	SIXT_ERROR("invalid XML loop parameters");
	return;
      }
      mydata->loop_depth++;
      return;

    } else {
      // Inner loop.
      mydata->further_loops=1;
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
  if (sprintf(buffer, "<%s", el)>=MAXMSG) {
    mydata->status=EXIT_FAILURE;
    SIXT_ERROR("XML element string too long");
    return;
  }
  addString2XMLBuffer(output, buffer, &mydata->status);
  CHECK_STATUS_VOID(mydata->status);

  int ii=0;
  while(attr[ii]) {
    if (sprintf(buffer, " %s=\"%s\"", attr[ii], attr[ii+1])>=MAXMSG) {
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
  int len_old=strlen(old);

  if ((0==strlen(old)) || (0==strlen(buffer->text))) return;

  while (NULL!=(occurence=strstr(buffer->text, old))) {
    // Get the length of the tail.
    int len_tail=strlen(occurence)-len_old;

    // String tail after the first occurence of the old string 
    // in the buffer text.
    char* tail=(char*)malloc((1+len_tail)*sizeof(char));
    if (NULL==tail) {
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for string buffer in XML "
		 "pre-parsing failed");
      return;
    }
    // Copy the tail of the string without the old string to the buffer.
    strcpy(tail, &occurence[len_old]);
    // Truncate the buffer string directly before the occurence 
    // of the old string.
    occurence[0]='\0';
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

  // Perform all "*" before "+" and "-" operations.
  // Loop while a "*" sign is found in the buffer text.
  while (NULL!=(occurrence=strpbrk(occurrence, "*"))) {

    // 1. Determine the first term.
    char* start=occurrence-1;
    // Scan forward until reaching a non-digit.
    while (strpbrk(start, "0123456789")==start) {
      start--;
    }
    start++;

    // Check if there is really a numeric term in front of the "+" or "-" sign.
    // If not continue with the next loop iteration.
    if (start==occurrence) {
      occurrence++;
      continue;
    }

    // Store the first value in a separate string.
    char svalue[MAXMSG];
    int ii;
    for (ii=0; start+ii<occurrence; ii++) {
      svalue[ii]=start[ii];
    }
    svalue[ii]='\0';
    
    // Convert the string to an integer value.
    int ivalue1=atoi(svalue);


    // 2. Determine the second term.
    char* end=occurrence+1;
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
      svalue[ii]=occurrence[1+ii];
    }
    svalue[ii]='\0';

    // Convert the string to an integer value.
    int ivalue2=atoi(svalue);


    // Perform the multiplication.
    int result=ivalue1*ivalue2;


    // Store the result at the right position in the XMLBuffer text.
    // Determine the new sub-string.
    sprintf(svalue, "%d", result);

    // Get the length of the tail.
    int len_tail=strlen(end)-1;

    // String tail after the first occurrence of the old string in the buffer text.
    char* tail=(char*)malloc((1+len_tail)*sizeof(char));
    if (NULL==tail) {
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for string buffer in XML pre-parser failed");
      return;
    }
    // Copy the tail of the string without the old string to the buffer.
    strcpy(tail, &end[1]);
    // Truncate the buffer string directly before the occurrence of the old string.
    start[0]='\0';
    // Append the new string to the truncated buffer string.
    addString2XMLBuffer(buffer, svalue, status);
    // Append the tail after the newly inserted new string.
    addString2XMLBuffer(buffer, tail, status);

    // Release memory of the tail buffer.
    free(tail);
  }
  // END of loop over all occurrences of "*" signs.

  // Reset the pointer to the beginning of the XML buffer text.
  occurrence=buffer->text;

  // Loop while a "+" or "-" sign is found in the buffer text.
  while (NULL!=(occurrence=strpbrk(occurrence, "+-"))) {

    // 1. Determine the first term.
    char* start=occurrence-1;
    // Scan forward until reaching a non-digit.
    while (strpbrk(start, "0123456789")==start) {
      start--;
    }
    start++;

    // Check if there is really a numeric term in front of the "+" or "-" sign.
    // If not (e.g. "e-4"), continue with the next loop iteration.
    if (start==occurrence) {
      occurrence++;
      continue;
    }

    // Store the first value in a separate string.
    char svalue[MAXMSG];
    int ii;
    for (ii=0; start+ii<occurrence; ii++) {
      svalue[ii]=start[ii];
    }
    svalue[ii]='\0';
    
    // Convert the string to an integer value.
    int ivalue1=atoi(svalue);


    // 2. Determine the second term.
    char* end=occurrence+1;
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
      svalue[ii]=occurrence[1+ii];
    }
    svalue[ii]='\0';

    // Convert the string to an integer value.
    int ivalue2=atoi(svalue);


    // Perform the arithmetic operation.
    int result=0;
    if (occurrence[0]=='+') {
      result=ivalue1+ivalue2;
    } else if (occurrence[0]=='-') {
      result=ivalue1-ivalue2;
    }


    // Store the result at the right position in the XMLBuffer text.
    // Determine the new sub-string.
    sprintf(svalue, "%d", result);

    // Get the length of the tail.
    int len_tail=strlen(end)-1;

    // String tail after the first occurrence of the old string 
    // in the buffer text.
    char* tail=(char*)malloc((1+len_tail)*sizeof(char));
    if (NULL==tail) {
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for string buffer in XML "
		 "pre-parsing failed");
      return;
    }
    // Copy the tail of the string without the old string to the buffer.
    strcpy(tail, &end[1]);
    // Truncate the buffer string directly before the occurrence 
    // of the old string.
    start[0]='\0';
    // Append the new string to the truncated buffer string.
    addString2XMLBuffer(buffer, svalue, status);
    // Append the tail after the newly inserted new string.
    addString2XMLBuffer(buffer, tail, status);

    // Release memory of the tail buffer.
    free(tail);
  }
  // END of loop over all occurrences of "+" or "-" signs.
}


static void expandXMLElementEnd(void* data, const char* el)
{
  struct XMLPreParseData* mydata=(struct XMLPreParseData*)data;

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
  if (sprintf(buffer, "</%s>", el)>=MAXMSG) {
    mydata->status=EXIT_FAILURE;
    SIXT_ERROR("XML string element too long");
    return;
  }
  addString2XMLBuffer(output, buffer, &mydata->status);
  if (EXIT_SUCCESS!=mydata->status) return;
}



static void InclXMLElementStart(void* data, 
			    const char* el, 
			    const char** attr) 
{
  struct XMLIncludeHandler* mydata=(struct XMLIncludeHandler*)data;

  // Check if an error has occurred previously.
  CHECK_STATUS_VOID(mydata->status);

  struct XMLBuffer* output=mydata->output_buffer;

  // Convert the element to an upper case string.
  char Uelement[MAXMSG]; // Upper case version of XML element
  strcpy(Uelement, el);
  strtoupper(Uelement);

  // Check if the element is an include tag.
  if(strcmp(Uelement, "INCLUDE")==0){

    // Set the output to the include buffer
    output=mydata->include_buffer;

    // Read the name of the included file
    char includefilename[MAXFILENAME], includefilepath[MAXFILENAME];
    getXMLAttributeString(attr, "FILENAME", includefilename);

    // Check if a file name has been specified.
    if(strlen(includefilename)==0){
      mydata->status=EXIT_FAILURE;
      SIXT_ERROR("Failed reading name of included xml file.");
      return;
    }

    // Construct the filepath to the included file
    strcpy(includefilepath, mydata->xmlfile);
    if(includefilename==NULL){
      mydata->status=EXIT_FAILURE;
      SIXT_ERROR("Failed copying name of included xml file.");
      return;
    }
    // Set a pointer to the last occurence of a slash
    char *ptr=strrchr(includefilepath, '/');
    if(ptr==NULL){
      // If no slash is found, set it to the first character
      ptr=&includefilepath[0];
    }else{
      // Otherwise move to the next character not to lose the slash
      ptr++;
    }
    // Set the next character to null to cut off the name of the input xml file
    *ptr='\0';
    // Append the name of the included xml file
    strcat(includefilepath, includefilename);

    // Open include file
    FILE* includefile=fopen(includefilepath, "r");
    if(includefile==NULL){
      mydata->status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "Failed opening included xml file:\n%s", includefilepath);
      SIXT_ERROR(msg);
      return;
    }
    // write include file to output buffer
    const int buffer_size=256;
    char buffer[buffer_size+1];
    int len;

    do{
      len=fread(buffer, 1, buffer_size, includefile);
      buffer[len]='\0';
      addString2XMLBuffer(output, buffer, &mydata->status);
      CHECK_STATUS_VOID(mydata->status);
    }while(!feof(includefile));

    fclose(includefile);

    // Recursively scan included xml code for includes
    expandIncludesXML(output, includefilepath, &mydata->status);
    CHECK_STATUS_VOID(mydata->status);

    // Copy included XML code to the right output buffer
    addString2XMLBuffer(mydata->output_buffer, output->text, &mydata->status);
    CHECK_STATUS_VOID(mydata->status);

    // raise further includes
    mydata->further_includes=1;

    // empty the include buffer
    freeXMLBuffer(&(mydata->include_buffer));
    mydata->include_buffer=newXMLBuffer(&mydata->status);

  }else{
    // If it was not an include tag, just write it to the output buffer
    char buffer[MAXMSG];
    if(sprintf(buffer, "<%s", el)>=MAXMSG){
      mydata->status=EXIT_FAILURE;
      SIXT_ERROR("XML element string too long");
      return;
    }
    addString2XMLBuffer(output, buffer, &mydata->status);
    CHECK_STATUS_VOID(mydata->status);
    int ii=0;
    while(attr[ii]){
      if(sprintf(buffer, " %s=\"%s\"", attr[ii], attr[ii+1])>=MAXMSG){
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
}

static void InclXMLElementEnd(void* data, 
			    const char* el) 
{
  struct XMLIncludeHandler* mydata=(struct XMLIncludeHandler*)data;

  // Pointer to output buffer
  struct XMLBuffer* output=mydata->output_buffer;

  // Convert the element to an upper case string.
  char Uelement[MAXMSG]; // Upper case version of XML element
  strcpy(Uelement, el);
  strtoupper(Uelement);

  // Check if the element is an include tag.
  if(strcmp(Uelement, "INCLUDE")){
    // If not, print the end tag
    char buffer[MAXMSG];
    if(sprintf(buffer, "</%s>", el)>=MAXMSG){
      mydata->status=EXIT_FAILURE;
      SIXT_ERROR("XML string element too long");
      return;
    }
    addString2XMLBuffer(output, buffer, &mydata->status);
    if(EXIT_SUCCESS!=mydata->status) {
      return;
    }
  }
}

void expandIncludesXML(struct XMLBuffer* const buffer, char* filename, int* const status)
{
  struct XMLIncludeHandler data;
  strcpy(data.xmlfile, filename);

  do{
    // Set further_includes to 0, if nothing new is found, 
    // this terminates the while loop
    data.further_includes=0;
    data.status=EXIT_SUCCESS;
    data.include_buffer=newXMLBuffer(status);
    data.output_buffer=newXMLBuffer(status);

    // Parse XML code in the buffer using the expat library.
    // Get a parser object.
    XML_Parser parser=XML_ParserCreate(NULL);
    if(NULL==parser){
      *status=EXIT_FAILURE;
      SIXT_ERROR("could not allocate memory for XML parser");
      return;
    }

    // Set data that is passed to the handler functions.
    XML_SetUserData(parser, &data);

    // Set handler functions
    XML_SetElementHandler(parser, InclXMLElementStart, InclXMLElementEnd);

    // Process all the data in the string buffer.
    const int done=1;
    if(!XML_Parse(parser, buffer->text, strlen(buffer->text), done)) {
      // Parse error.
      *status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "parsing XML code failed: \n%s\n", 
	      XML_ErrorString(XML_GetErrorCode(parser)));
      printf("%s", buffer->text);
      SIXT_ERROR(msg);
      return;
    }

    // Check for errors.
    if (EXIT_SUCCESS!=data.status) {
      *status=data.status;
      return;
    }

    // Copy the output XMLBuffer to the input XMLBuffer
    copyXMLBuffer(buffer, data.output_buffer, status);
    if(EXIT_SUCCESS!=*status){
      return;
    }
    // release allocated memory
    freeXMLBuffer(&data.output_buffer);
    XML_ParserFree(parser);

  }while(data.further_includes);

  CHECK_STATUS_VOID(*status);
}


void expandXML(struct XMLBuffer* const buffer, int* const status)
{
  struct XMLPreParseData data;
  
  do {

    // Parse XML code in the xmlbuffer using the expat library.
    // Get an XML_Parser object.
    XML_Parser parser=XML_ParserCreate(NULL);
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
    data.further_loops =0;
    data.loop_depth    =0;
    data.loop_start    =0;
    data.loop_end      =0;
    data.loop_increment=0;
    data.output_buffer =newXMLBuffer(status);
    data.loop_buffer   =newXMLBuffer(status);
    data.status=EXIT_SUCCESS;

    // Process all the data in the string buffer.
    const int done=1;
    if (!XML_Parse(parser, buffer->text, strlen(buffer->text), done)) {
      // Parse error.
      *status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "parsing XML code failed:\n%s\n", 
	      XML_ErrorString(XML_GetErrorCode(parser)));
      printf("%s", buffer->text);
      SIXT_ERROR(msg);
      return;
    }
    // Check for errors.
    if (EXIT_SUCCESS!=data.status) {
      *status=data.status;
      return;
    }

    // Copy the output XMLBuffer to the input XMLBuffer ...
    copyXMLBuffer(buffer, data.output_buffer, status);
    if (EXIT_SUCCESS!=*status) return;
    // ... and release allocated memory.
    freeXMLBuffer(&data.output_buffer);
    freeXMLBuffer(&data.loop_buffer);

    XML_ParserFree(parser);

  } while (data.further_loops>0);

  // Replace arithmetic +/- expressions.
  execArithmeticOpsInXMLBuffer(buffer, status);
  CHECK_STATUS_VOID(*status);

  // Replace all '\', which are used to escape e.g. '-' signs in
  // file names in order to avoid misplaced arithmetic evaluations.
  replaceInXMLBuffer(buffer, "\\", "", status);
  CHECK_STATUS_VOID(*status);
}


void getXMLAttributeString(const char** attr, 
			   const char* const key, 
			   char* const value)
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


float getXMLAttributeFloat(const char** attr, const char* const key)
{
  char buffer[MAXMSG]; // String buffer.
  getXMLAttributeString(attr, key, buffer);
  return((float)atof(buffer));
}


double getXMLAttributeDouble(const char** attr, const char* const key)
{
  char buffer[MAXMSG]; // String buffer.
  getXMLAttributeString(attr, key, buffer);
  return(atof(buffer));
}


int getXMLAttributeInt(const char** attr, const char* const key)
{
  char buffer[MAXMSG]; // String buffer.
  getXMLAttributeString(attr, key, buffer);
  return(atoi(buffer));
}


long getXMLAttributeLong(const char** attr, const char* const key)
{
  char buffer[MAXMSG]; // String buffer.
  getXMLAttributeString(attr, key, buffer);
  return(atol(buffer));
}


