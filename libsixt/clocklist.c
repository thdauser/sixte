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

  // Set the initial values.
  list->nelements=0;
  list->element  =0;
  list->time     =0.;
  list->readout_time=0.;
  list->frame    =0;

  return(list);
}


void destroyClockList(ClockList** const list)
{
  if (NULL!=(*list)) {
    if ((NULL!=(*list)->list)&&(NULL!=(*list)->type)&&(0<(*list)->nelements)) {
      int i;
      for (i=0; i<(*list)->nelements; i++) {
	if (NULL!=(*list)->list[i]) {
	  if (CL_WAIT==(*list)->type[i]) {
	    CLWait* clwait=(CLWait*)((*list)->list[i]);
	    destroyCLWait(&clwait);
	    (*list)->list[i]=NULL;
	  } else if (CL_LINESHIFT==(*list)->type[i]) {
	    CLLineShift* cllineshift = (CLLineShift*)((*list)->list[i]);
	    destroyCLLineShift(&cllineshift);
	    (*list)->list[i]=NULL;
	  } else if (CL_NEWFRAME==(*list)->type[i]) {
	    CLNewFrame* clnewframe=(CLNewFrame*)((*list)->list[i]);
	    destroyCLNewFrame(&clnewframe);
	    (*list)->list[i]=NULL;
	  } else if (CL_READOUTLINE==(*list)->type[i]) {
	    CLReadoutLine* clreadoutline=(CLReadoutLine*)((*list)->list[i]);
	    destroyCLReadoutLine(&clreadoutline);
	    (*list)->list[i]=NULL;
	  } else if (CL_CLEARLINE==(*list)->type[i]) {
	    CLClearLine* clclearline=(CLClearLine*)((*list)->list[i]);
	    destroyCLClearLine(&clclearline);
	    (*list)->list[i]=NULL;
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
  list->type=realloc(list->type, (list->nelements+1)*sizeof(CLType));
  if (NULL==list->type) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for ClockList element failed");
    return;
  }
  list->list = realloc(list->list, (list->nelements+1)*sizeof(void*));
  if (NULL==list->list) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for ClockList element failed");
    return;
  }
  list->nelements++;

  // Setup the values of the last, newly appended element.
  list->type[list->nelements-1]=type;
  list->list[list->nelements-1]=element;
}



static inline void moveClockList2NextElement(ClockList* const list)
{
  list->element++;

  // Check if the end of the list is reached.
  if (list->element>=list->nelements) {
    // Start from the beginning.
    list->element=0;
  }
}


void getClockListElement(ClockList* const list, const double time,
			 CLType* const type, void** const element,
			 int* const status)
{
  // Check if the list contains any elements.
  if (0==list->nelements) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("clock list is empty");
    return;
  }

  // If the current element is of the type CL_WAIT, we have to check
  // whether the time specified in the parameter list exceeds the waiting
  // period. In that case move to the next element.
  if (CL_WAIT==list->type[list->element]) {
    CLWait* clwait=(CLWait*)list->list[list->element];
    if (list->time+clwait->time>=time) {
      // The wait period still continues.
      *type   =CL_NONE;
      *element=NULL;
      return;
    } else {
      // The wait period is finished.
      list->time+=clwait->time;
    }

  } else if (CL_NEWFRAME==list->type[list->element]) {
    // Start a new frame.
    list->frame++;
    list->readout_time=list->time;
  }

  *type   =list->type[list->element];
  *element=list->list[list->element];
  moveClockList2NextElement(list);
}


CLWait* newCLWait(const double time, int* const status) {
  headas_chat(5, "new CLWait element (time=%e)\n", time);

  // Allocate memory.
  CLWait* clwait=(CLWait*)malloc(sizeof(CLWait));
  if (NULL==clwait) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for CLWait element failed");
    return(clwait);
  }

  // Initialize pointers with NULL.

  // Set parameters.
  clwait->time=time;

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
  CLLineShift* cllineshift=(CLLineShift*)malloc(sizeof(CLLineShift));
  if (NULL==cllineshift) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for CLLineShift element failed");
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


CLNewFrame* newCLNewFrame(int* const status) {
  headas_chat(5, "new CLNewFrame element\n");

  // Allocate memory.
  CLNewFrame* clnewframe=(CLNewFrame*)malloc(sizeof(CLNewFrame));
  if (NULL==clnewframe) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for CLNewFrame element failed");
    return(clnewframe);
  }

  return(clnewframe);
}


void destroyCLNewFrame(CLNewFrame** const clnewframe)
{
  if (NULL!=(*clnewframe)) {
    free(*clnewframe);
    *clnewframe=NULL;
  }
}


CLReadoutLine* newCLReadoutLine(const int lineindex, const int readoutindex,
				int* const status)
{
  headas_chat(5, "new CLReadoutLine element (lineindex=%d, readoutindex=%d)\n",
	      lineindex, readoutindex);

  // Allocate memory.
  CLReadoutLine* clreadoutline=(CLReadoutLine*)malloc(sizeof(CLReadoutLine));
  if (NULL==clreadoutline) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for CLReadoutLine element failed");
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


CLClearLine* newCLClearLine(const int lineindex, int* const status)
{
  headas_chat(5, "new CLClearLine element (lineindex=%d)\n",
	      lineindex);

  // Allocate memory.
  CLClearLine* clclearline=(CLClearLine*)malloc(sizeof(CLClearLine));
  if (NULL==clclearline) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for CLClearLine element failed");
    return(clclearline);
  }

  // Initialize.
  clclearline->lineindex=lineindex;

  return(clclearline);
}



void destroyCLClearLine(CLClearLine** const clclearline)
{
  if (NULL!=(*clclearline)) {
    free(*clclearline);
    *clclearline=NULL;
  }
}
