#include "xmlbuffer.h"


void addString2XMLBuffer(struct XMLBuffer* const buffer, 
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


struct XMLBuffer* newXMLBuffer(int* const status)
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

  // Perform all "*" before "+" and "-" operations.
  // Loop while a "*" sign is found in the buffer text.
  while (NULL!=(occurrence=strpbrk(occurrence, "*"))) {

    // 1. Determine the first term.
    char* start = occurrence-1;
    // Scan forward until reaching a non-digit.
    while (strpbrk(start, "0123456789")==start) {
      start--;
    }
    start++;

    // Check if there really is a numeric term in front of the "+" or "-" sign.
    // If not continue with the next loop iteration.
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


    // Perform the multiplication.
    int result = ivalue1 * ivalue2;


    // Store the result at the right position in the XMLBuffer text.
    // Determine the new sub-string.
    sprintf(svalue, "%d", result);

    // Get the length of the tail.
    int len_tail = strlen(end)-1;

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
    start[0] = '\0';
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

	freeXMLBuffer(&replacedBuffer);
      }
      // Clear the loop buffer.
      freeXMLBuffer(&mydata->loop_buffer);
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


void expandXML(struct XMLBuffer* const buffer, int* const status)
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
    freeXMLBuffer(&data.output_buffer);
    freeXMLBuffer(&data.loop_buffer);

    XML_ParserFree(parser);

  } while (data.further_loops>0);

  // Replace arithmetic +/- expressions.
  execArithmeticOpsInXMLBuffer(buffer, status);
}
