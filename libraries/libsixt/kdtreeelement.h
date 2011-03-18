#ifndef KDTREEELEMENT_H
#define KDTREEELEMENT_H 1

#include "sixt.h"
#include "source.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Element of a KDTree (multidimensional binary tree). */
struct structKDTreeElement {
  Source* src;
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

/** Build up the KDTree from the given list of Sources. */
KDTreeElement* buildKDTree2(Source* const list, 
			    const long nelements,
			    const int depth,
			    int* const status);

/** Perform a range search on the given kdTree, i.e., return all X-ray
    sources lying within a certain radius around the reference
    point. This region is defined by the minimum cosine value for the
    scalar product of the source direction and the reference
    vector. The function returns a time-ordered list of newly
    generated photons. */
LinkedPhoListElement* KDTreeRangeSearch(KDTreeElement* const node, 
					const int depth,
					const Vector* const ref, 
					const double min_align, 
					const double t0, const double t1,
					int* const status);


#endif /* KDTREEELEMENT_H */
