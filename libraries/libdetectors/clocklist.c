#include "clocklist.h"


ClockList* newClockList(int* const status)
{
  headas_chat(5, "initialize empty ClockList object ...\n");

  // Allocate memory.
  ClockList* list = (ClockList*)malloc(sizeof(ClockList));
  if (NULL==list) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for ClockList failed!\n", *status);
    return(list);
  }

  // Initialize all pointers with NULL.
  list->type=NULL;
  list->list=NULL;
  list->nelements=0;
  list->element  =0;
  list->time     =0.;

  return(list);
}



void destroyClockList(ClockList** const list)
{
  if (NULL!=(*list)) {
    if ( (NULL!=(*list)->list) && (NULL!=(*list)->type) && (0<(*list)->nelements) ) {
      int i;
      for (i=0; i<(*list)->nelements; i++) {
	if (NULL!=(*list)->list[i]) {
	  if (CL_WAIT==(*list)->type[i]) {
	    CLWait* clwait = (CLWait*)((*list)->list[i]);
	    destroyCLWait(&clwait);
	    (*list)->list[i] = NULL;
	  } else if (CL_LINESHIFT==(*list)->type[i]) {
	    CLLineShift* cllineshift = (CLLineShift*)((*list)->list[i]);
	    destroyCLLineShift(&cllineshift);
	    (*list)->list[i] = NULL;
	  } else if (CL_READOUTLINE==(*list)->type[i]) {
	    CLReadoutLine* clreadoutline = (CLReadoutLine*)((*list)->list[i]);
	    destroyCLReadoutLine(&clreadoutline);
	    (*list)->list[i] = NULL;
	  } else {
	    printf("Unknown CLType!\n");
	  }
	}
      }
      // END of loop over all elements in the list.
      free((*list)->list);
      free((*list)->type);
    }
    free(*list);
    *list=NULL;
  }
}



void append2ClockList(ClockList* const list, const CLType type, 
		      void* const element, int* const status)
{
  // Allocate memory for a new element.
  list->type = realloc(list->type, (list->nelements+1)*sizeof(CLType));
  if (NULL==list->type) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for ClockList element failed!\n", *status);
    return;
  }
  list->list = realloc(list->list, (list->nelements+1)*sizeof(void*));
  if (NULL==list->list) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for ClockList element failed!\n", *status);
    return;
  }
  list->nelements++;
  
  // Setup the values of the last, newly appended element.
  list->type[list->nelements-1] = type;
  list->list[list->nelements-1] = element;
}



static inline void moveClockList2NextElement(ClockList* const list)
{
  list->element++;
  if (list->element>=list->nelements) {
    list->element=0;
  }
}



void getClockListElement(ClockList* const list, const double time,
			 CLType* type, void** const element)
{
  // Check if the list contains any elements.
  if (0==list->nelements) {
    type=CL_NONE;
    return;
  }

  // If the current element is of the type CL_WAIT, we have to check
  // whether the time specified in the parameter list exceeds the waiting
  // period. In that case move to the next element.
  while (CL_WAIT==list->type[list->element]) {
    CLWait* clwait = (CLWait*)list->list[list->element];
    if (list->time+clwait->time >= time) {
      *type   =CL_WAIT;
      *element=list->list[list->element];
      return;
    } else {
      moveClockList2NextElement(list);
      list->time+=clwait->time;
    }
  }
  
  // The current element is NOT of type CL_WAIT.
  *type   =list->type[list->element];
  *element=list->list[list->element];
  moveClockList2NextElement(list);
}



CLWait* newCLWait(const double time, int* const status) {
  headas_chat(5, "new CLWait element (time=%e)\n", time);

  // Allocate memory.
  CLWait* clwait = (CLWait*)malloc(sizeof(CLWait));
  if (NULL==clwait) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for CLWait element failed!\n", *status);
    return(clwait);
  }
  
  // Initialize.
  clwait->time = time;

  return(clwait);
}



void destroyCLWait(CLWait** const clwait)
{
  if (NULL!=(*clwait)) {
    free(*clwait);
    *clwait=NULL;
  }
}



CLLineShift* newCLLineShift(int* const status) {
  headas_chat(5, "new CLLineShift element\n");

  // Allocate memory.
  CLLineShift* cllineshift = (CLLineShift*)malloc(sizeof(CLLineShift));
  if (NULL==cllineshift) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for CLLineShift element failed!\n", *status);
    return(cllineshift);
  }
  
  return(cllineshift);
}



void destroyCLLineShift(CLLineShift** const cllineshift)
{
  if (NULL!=(*cllineshift)) {
    free(*cllineshift);
    *cllineshift=NULL;
  }
}



CLReadoutLine* newCLReadoutLine(const int lineindex, const int readoutindex, 
				int* const status)
{
  headas_chat(5, "new CLReadoutLine element (lineindex=%d, readoutindex=%d)\n", 
	      lineindex, readoutindex);

  // Allocate memory.
  CLReadoutLine* clreadoutline = (CLReadoutLine*)malloc(sizeof(CLReadoutLine));
  if (NULL==clreadoutline) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for CLReadoutLine element failed!\n", *status);
    return(clreadoutline);
  }
  
  // Initialize.
  clreadoutline->lineindex    = lineindex;
  clreadoutline->readoutindex = readoutindex;

  return(clreadoutline);
}



void destroyCLReadoutLine(CLReadoutLine** const clreadoutline)
{
  if (NULL!=(*clreadoutline)) {
    free(*clreadoutline);
    *clreadoutline=NULL;
  }
}

