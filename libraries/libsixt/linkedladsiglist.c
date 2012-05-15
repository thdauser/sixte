#include "linkedladsiglist.h"

LinkedLADSigListElement* newLinkedLADSigListElement(int* const status)
{
  LinkedLADSigListElement* el= 
    (LinkedLADSigListElement*)malloc(sizeof(LinkedLADSigListElement));
  CHECK_NULL(el, *status,
	     "memory allocation for LinkedLADSigListElement failed");
  
  // Initialize pointers with NULL.
  el->next=NULL;

  return(el);
}


void freeLinkedLADSigList(LinkedLADSigListElement** const list)
{
  if (NULL!=*list) {
    if (NULL!=(*list)->next) {
      freeLinkedLADSigList(&((*list)->next));
    }

    free(*list);
    *list=NULL;
  }
}
