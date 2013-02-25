#include "linkedpholist.h"

LinkedPhoListElement* newLinkedPhoListElement(int* const status)
{
  LinkedPhoListElement* el=
    (LinkedPhoListElement*)malloc(sizeof(LinkedPhoListElement));
  CHECK_NULL(el, *status,
	     "memory allocation for LinkedPhoListElement failed");
  
  // Initialize pointers with NULL.
  el->next=NULL;

  return(el);
}


void freeLinkedPhoList(LinkedPhoListElement** const list)
{
  if (NULL!=*list) {
    if (NULL!=(*list)->next) {
      freeLinkedPhoList(&((*list)->next));
    }

    free(*list);
    *list=NULL;
  }
}


LinkedPhoListElement* mergeLinkedPhoLists(LinkedPhoListElement* list1, 
					  LinkedPhoListElement* list2)
{
  LinkedPhoListElement* start   =NULL;
  LinkedPhoListElement** current=&start;

  if (NULL==list1) return(list2);
  if (NULL==list2) return(list1);

  while((NULL!=list1) && (NULL!=list2)) {
    if (list1->photon.time <= list2->photon.time) {
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

