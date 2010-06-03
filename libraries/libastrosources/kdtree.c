#include "kdtree.h"


kdTree buildKDTree(PointSource* list, long nelements)
{
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
		       Vector* ref, double radius2, 
		       double time, double dt, 
		       struct PhotonOrderedListEntry** pl,
		       struct RMF* rmf,
		       int* status)
{
  // Check if the kd-Tree exists.
  if (NULL==node) return;

  // Calculate the distance (squared) between the node and the 
  // reference point.
  Vector location = unit_vector(node->source.ra, node->source.dec);
  double distance2 = 
    pow(location.x-ref->x, 2.) +
    pow(location.y-ref->y, 2.) +
    pow(location.z-ref->z, 2.);

  // Check if the current node lies within the search radius.
  if (distance2 <= radius2) {
    create_PointSourcePhotons(&node->source, time, dt, pl, rmf);
  }

  // Check if we are at a leaf.
  if ((NULL==node->left) && (NULL==node->right)) return;

  int axis = depth % 3;

  // Check which branch to search first.
  kdNode* near;
  kdNode* far;
  double distance2edge = 
    getVectorDimensionValue(ref, axis) -
    getVectorDimensionValue(&location, axis);
  if (distance2edge < 0.) {
    near = node->left;
    far  = node->right;
  } else {
    far  = node->left;
    near = node->right;
  }

  // Descent into near tree if it exists, and then check
  // against current node.
  if (NULL!=near) {
    kdTreeRangeSearch(near, depth+1, ref, radius2,
		      time, dt, pl, rmf, status);
  } 
  // END of (NULL!=near)

  // Check whether we have to look into the far tree.
  // A search is only necessary if the minimum distance
  // of the reference point is such that we can have an
  // overlap there.
  if (NULL!=far) {
    if (distance2edge*distance2edge < radius2) {
      // Move to the end of the linked list.
      // Append newly found entries.
      kdTreeRangeSearch(far, depth+1, ref, radius2,
			time, dt, pl, rmf, status);
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
    free(tree);
  }
}

