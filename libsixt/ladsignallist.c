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
