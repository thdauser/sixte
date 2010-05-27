#include "kdtree.h"


kdNode* kdTreeBuild(SourceList* list, long nelements, int depth)
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
    node->location = list[0].location;
    node->left  = NULL;
    node->right = NULL;
    return(node);
  }

  long median = nelements/2;
  int axis = depth % 3;
  quicksortSourceList(list, 0, nelements-1, axis);

  // Fill the newly created node with data.
  node->location = list[median].location;

  // Set right and left pointers of node.
  if (median>0) {
    node->left = kdTreeBuild(list, median, depth+1);
  } else {
    node->left = NULL;
  }

  if (median<nelements-1) {
    node->right = kdTreeBuild(&list[median+1], 
			      nelements-median-1, 
			      depth+1);
  } else {
    node->right = NULL;
  }

  return(node);
}



int addNode2SourceList(kdNode* node, SourceList** list, long* nelements)
{
  // Check if new memory has to be allocated.
  (*nelements)++;
  if (1 == (*nelements % 1000)) {
    *list = (SourceList*)realloc(*list, ((*nelements/1000)+1)*sizeof(SourceList));
    if (NULL==*list) {
      HD_ERROR_THROW("Error: Could not allocate memory for SourceList!\n",
		     EXIT_FAILURE);
      return(EXIT_FAILURE);
    }
  }

  // Add the new Source to the list.
  (*list)[(*nelements)-1].location = node->location;

  return(EXIT_SUCCESS);
}



int kdTreeRangeSearch(kdNode* node, int depth,
		      Vector* ref, double radius2, 
		      SourceList** list, long *nelements)
{
  int status = EXIT_SUCCESS;

  // Calculate the distance (squared) between the node and the 
  // reference point.
  double distance2 = 
    pow(node->location.x-ref->x, 2.) +
    pow(node->location.y-ref->y, 2.) +
    pow(node->location.z-ref->z, 2.);

  // Check if the current node lies within the search radius.
  if (distance2 <= radius2) {
    status = addNode2SourceList(node, list, nelements);
    if (EXIT_SUCCESS!=status) return(status);
  }

  // Check if we are at a leaf.
  if ((NULL==node->left) && (NULL==node->right)) return(status);

  int axis = depth % 3;

  // Check which branch to search first.
  kdNode* near;
  kdNode* far;
  double distance2edge = 
    getVectorDimensionValue(ref, axis) -
    getVectorDimensionValue(&node->location, axis);
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
		      list, nelements);
  } 
  // END of (NULL!=near)

  // Check whether we have to look into the far tree.
  // A search is only necessary if the minimum distance
  // of the reference point is such that we can have an
  // overlap there.
  if (NULL!=far) {
    if (distance2edge*distance2edge < radius2) {
      kdTreeRangeSearch(far, depth+1, ref, radius2,
			list, nelements);
    }
  }
  // END of (NULL!=far)

  return(status);
}




