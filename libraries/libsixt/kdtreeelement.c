#include "kdtreeelement.h"


KDTreeElement* newKDTreeElement(int* const status)
{
  KDTreeElement* el = (KDTreeElement*)malloc(sizeof(KDTreeElement));
  CHECK_NULL(el, *status, 
	     "memory allocation for KDTreeElement failed");

  // Initalize pointers with NULL.
  el->element=NULL;
  el->left   = NULL;
  el->right  = NULL;

  return(el);
}


void freeKDTreeElement(KDTreeElement** el)
{
  if (NULL!=*el) {
    if (NULL!=(*el)->left) {
      freeKDTreeElement(&((*el)->left));
    }
    if (NULL!=(*el)->right) {
      freeKDTreeElement(&((*el)->right));
    }
    assert(NULL==(*el)->element);

    free(*el);
    *el=NULL;
  }
}

