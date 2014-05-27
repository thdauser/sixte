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
*/

#include "ladsignallist.h"

LADSignalListItem* newLADSignalListItem(int* const status)
{
  LADSignalListItem* el= 
    (LADSignalListItem*)malloc(sizeof(LADSignalListItem));
  CHECK_NULL(el, *status,
	     "memory allocation for LADSignalListItem failed");
  
  // Initialize pointers with NULL.
  el->next=NULL;

  return(el);
}


void freeLADSignalList(LADSignalListItem** const list)
{
  if (NULL!=*list) {
    if (NULL!=(*list)->next) {
      freeLADSignalList(&((*list)->next));
    }

    free(*list);
    *list=NULL;
  }
}
