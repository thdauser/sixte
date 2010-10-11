#include "kdtree.h"


kdTree buildKDTree(PointSource* list, long nelements)
{
  //printf("Build KDTree with %ld elements\n", nelements);
  kdnelements=nelements;

  // Check if the list is empty.
  if (0==nelements) return(NULL);

  return(buildKDNode(list, nelements, 0));
}



kdNode* buildKDNode(PointSource* list, long nelements, int depth)
{
  // Get a new empty node.
  kdNode* node = (kdNode*)malloc(sizeof(kdNode));
  if (NULL==node) {
    HD_ERROR_THROW("Error: Could not allocate memory for kdNode!", 
		   EXIT_FAILURE);
    return(node);
  };

  // Check if there is only one element in the source list.
  if (1==nelements) {
    node->source = list[0];
    node->left   = NULL;
    node->right  = NULL;
    return(node);
  }

  long median = nelements/2;
  int axis = depth % 3;
  quicksortPointSources(list, 0, nelements-1, axis);

  // Fill the newly created node with data.
  node->source = list[median];

  // Set right and left pointers of node.
  if (median>0) {
    node->left = buildKDNode(list, median, depth+1);
  } else {
    node->left = NULL;
  }

  if (median<nelements-1) {
    node->right = buildKDNode(&list[median+1], 
			      nelements-median-1, 
			      depth+1);
  } else {
    node->right = NULL;
  }

  return(node);
}



void kdTreeRangeSearch(kdNode* node, int depth,
		       Vector* ref, double min_align, 
		       double time, double dt, 
		       struct PhotonOrderedListEntry** pl,
		       const struct ARF* const arf,
		       int* status)
{
  // Check if the kd-Tree exists.
  if (NULL==node) return;

  // Check if the current node lies within the search radius.
  kdnchecked++; // RM
  if (fabs(scalar_product(&node->source.location, ref)) > min_align) {
    kdnfound++;
    create_PointSourcePhotons(&node->source, time, dt, pl, arf);
  }

  // Check if we are at a leaf.
  if ((NULL==node->left) && (NULL==node->right)) {
    return;
  }

  int axis = depth % 3;

  // Check which branch to search first.
  kdNode* near;
  kdNode* far;
  double distance2edge = 
    getVectorDimensionValue(ref, axis) -
    getVectorDimensionValue(&node->source.location, axis);
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
    kdTreeRangeSearch(near, depth+1, ref, min_align,
		      time, dt, pl, arf, status);
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
      kdTreeRangeSearch(far, depth+1, ref, min_align,
			time, dt, pl, arf, status);
    }
  }
  // END of (NULL!=far)

  return;
}



void freeKDTree(kdNode* tree)
{
  if (NULL!=tree) {
    freeKDTree(tree->left);
    freeKDTree(tree->right);
    freePointSource(tree->source);
    free(tree);
  }
}

