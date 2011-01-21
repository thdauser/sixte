#ifndef LINKEDLIST_H
#define LINKEDLIST_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Linked list. */
struct structLinkedListElement {

  void* el;

  /** Next entry. */
  struct structLinkedListElement* next;
};
typedef struct structLinkedListElement LinkedListElement;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. */
LinkedListElement* newLinkedListElement(int* const status);

/** Destructor. */
void freeLinkedList(LinkedListElement** list);

/** Append a new element at the end of the linked list. The function
    returns a pointer to the newly created last element in the linked
    list. */
LinkedListElement* appendEl2LinkedList(LinkedListElement** list, void* el, 
				       int* const status);


#endif /* LINKEDLIST_H */
