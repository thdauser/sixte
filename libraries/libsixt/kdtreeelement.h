#ifndef KDTREEELEMENT_H
#define KDTREEELEMENT_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Element of a KDTree (multidimensional binary tree). */
struct structKDTreeElement {
  void* element;
  struct structKDTreeElement* left;
  struct structKDTreeElement* right;
};
typedef struct structKDTreeElement KDTreeElement;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. */
KDTreeElement* newKDTreeElement(int* const status);

/** Destructor. */
void freeKDTreeElement(KDTreeElement** el);


#endif /* KDTREEELEMENT_H */
