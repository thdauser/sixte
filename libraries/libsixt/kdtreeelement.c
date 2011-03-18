#include "kdtreeelement.h"


KDTreeElement* newKDTreeElement(int* const status)
{
  KDTreeElement* el = (KDTreeElement*)malloc(sizeof(KDTreeElement));
  CHECK_NULL(el, *status, 
	     "memory allocation for KDTreeElement failed");

  // Initalize pointers with NULL.
  el->src   = NULL;
  el->left  = NULL;
  el->right = NULL;

  // Get memory for the content.
  el->src = newSource(status);
  CHECK_STATUS_RET(*status, el);

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
    freeSource(&(*el)->src);

    free(*el);
    *el=NULL;
  }
}


KDTreeElement* buildKDTree2(Source* const list, 
			    const long nelements,
			    const int depth,
			    int* const status)
{
  if (0==nelements) return(NULL);

  // Get a new empty node.
  KDTreeElement* node = newKDTreeElement(status);
  CHECK_STATUS_RET(*status, node);

  // Check if there is only one element in the source list.
  if (1==nelements) {
    *(node->src) = list[0];
    return(node);
  }

  long median = nelements/2;
  int axis = depth % 3;
  quicksortSources(list, 0, nelements-1, axis);

  // Fill the newly created node with data.
  *(node->src) = list[median];

  // Set right and left pointers of node.
  if (median>0) {
    node->left = buildKDTree2(list, median, depth+1, status);
  }

  if (median<nelements-1) {
    node->right = buildKDTree2(&list[median+1], 
			       nelements-median-1, 
			       depth+1, status);
  }

  return(node);
}



LinkedPhoListElement* KDTreeRangeSearch(KDTreeElement* const node, 
					const int depth,
					const Vector* const ref, 
					const double min_align, 
					const double t0, const double t1,
					int* const status)
{
  // Check if the kd-Tree exists.
  if (NULL==node) return(NULL);
  
  LinkedPhoListElement* list = NULL;
  
  // Check if the current node lies within the search radius.
  Vector location = unit_vector(node->src->ra, node->src->dec);
  if (fabs(scalar_product(&location, ref)) > min_align) {
    // Generate photons for this particular source.
    list = getXRayPhotons(node->src, t0, t1, status);
    CHECK_STATUS_RET(*status, list);
  }

  // Check if we are at a leaf.
  if ((NULL==node->left) && (NULL==node->right)) {
    return(list);
  }

  int axis = depth % 3;

  // Check which branch to search first.
  KDTreeElement* near;
  KDTreeElement* far;
  double distance2edge = 
    getVectorDimensionValue(ref, axis)-getVectorDimensionValue(&location, axis);
  if (distance2edge < 0.) {
    near = node->left;
    far  = node->right;
  } else {
    far  = node->left;
    near = node->right;
  }

  // Descend into near tree if it exists, and then check
  // against current node.
  if (NULL!=near) {
    LinkedPhoListElement* near_list = 
      KDTreeRangeSearch(near, depth+1, ref, min_align, t0, t1, status);

    // Merge the new photons into the existing list.
    list = mergeLinkedPhoLists(list, near_list);
  } 
  // END of (NULL!=near)

  // Check whether we have to look into the far tree.
  // A search is only necessary if the minimum distance
  // of the reference point is such that we can have an
  // overlap there.
  if (NULL!=far) {
    if (cos(distance2edge) > min_align) {
      // Move to the end of the linked list.
      // Append newly found entries.
      LinkedPhoListElement* far_list = 
	KDTreeRangeSearch(far, depth+1, ref, min_align, t0, t1, status);

      // Merge the new photons into the existing list.
      list = mergeLinkedPhoLists(list, far_list);
    }
  }
  // END of (NULL!=far)

  return(list);
}

