#include "linkedlist.h"

LinkedListElement* newLinkedListElement(int* const status)
{
  LinkedListElement* el = (LinkedListElement*)malloc(sizeof(LinkedListElement));
  CHECK_NULL(el, *status,
	     "memory allocation for LinkedListElement failed");
  
  // Initialize pointers with NULL.
  el->el   = NULL;
  el->next = NULL;

  return(el);
}


void freeLinkedList(LinkedListElement** list)
{
  if (NULL!=*list) {
    if (NULL!=(*list)->next) {
      freeLinkedList(&((*list)->next));
    }
    assert(NULL==(*list)->el);

    free(*list);
    *list=NULL;
  }
}


