#ifndef KDTREE_H
#define KDTREE_H 1

#include "sixt.h"
#include "vector.h"
#include "pointsources.h"
#include "pointsourcelist.h"
#include "photon.h"


#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////

long kdnelements;// Zahl der Elemente im KDTree.
long kdnfound;   // Zahl der gefundenen Einträge.
long kdnchecked; // Zahl der überprüften Einträge.

/** Node in the kdTree containing the X-ray sources in 3-dimensional
    space. */
struct structkdNode {
  /** Data structure representing the X-ray source. */
  PointSource source;

  struct structkdNode* left; /**< Pointer to node on the left. */
  struct structkdNode* right; /**< Pointer to node on the right. */
};
typedef struct structkdNode kdNode;

typedef kdNode* kdTree;


////////////////////////////////////////////////////////////////////////
//   Function declarations
////////////////////////////////////////////////////////////////////////


/** Build up the kdTree from the given list of sources. */
kdTree buildKDTree(PointSource* list, long nelements);

kdNode* buildKDNode(PointSource* list, long nelements, int depth);

/** Perform a range search on the given kbTree, i.e., return all X-ray
    sources lying within a certain radius around the reference
    point. This region is defined by the minimum cosine value for the
    scalar product of the source direction and the reference
    vector. New sources are appended at the end of the SourceList and
    the number of the returned X-ray sources is stored in the
    nelements parameter. The function return value is the error
    status. */
void kdTreeRangeSearch(kdNode* node, int depth,
		       Vector* ref, double min_align, 
		       double time, double dt, 
		       struct PhotonOrderedListEntry** list_first,
		       struct RMF* rmf,
		       int* status);

/** Destructor. */
void freeKDTree(kdNode* tree);


#endif /* KDTREE_H */
