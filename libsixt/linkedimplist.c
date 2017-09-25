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

#include "linkedimplist.h"

LinkedImpListElement* newLinkedImpListElement(int* const status){
	LinkedImpListElement* el=(LinkedImpListElement*)malloc(sizeof(LinkedImpListElement));
	CHECK_NULL(el, *status,
	     "memory allocation for LinkedImpListElement failed");

	// Initialize pointers with NULL.
	el->next=NULL;

	return(el);
}

void freeLinkedImpList(LinkedImpListElement** const list)
{
  if (NULL!=*list) {
    if (NULL!=(*list)->next) {
      freeLinkedImpList(&((*list)->next));
    }

    free(*list);
    *list=NULL;
  }
}


LinkedImpListElement* mergeLinkedImpLists(LinkedImpListElement* list1,
					  LinkedImpListElement* list2)
{
  LinkedImpListElement* start   =NULL;
  LinkedImpListElement** current=&start;

  if (NULL==list1) return(list2);
  if (NULL==list2) return(list1);

  while((NULL!=list1) && (NULL!=list2)) {
    if (list1->impact.time <= list2->impact.time) {
      *current = list1;
      list1 = list1->next;
    } else {
      *current = list2;
      list2 = list2->next;
    }
    current = &((*current)->next);
  }
  if (NULL!=list1) {
    *current = list1;
  }
  if (NULL!=list2) {
    *current = list2;
  }

  return(start);
}

