#include "sourcecatalog.h"


SourceCatalog* newSourceCatalog(int* const status)
{
  SourceCatalog* cat = (SourceCatalog*)malloc(sizeof(SourceCatalog));
  CHECK_NULL(cat, *status,
	     "memory allocation for SourceCatalog failed");

  // Initialize pointers with NULL.
  cat->tree   = NULL;
  cat->simput = NULL;

  return(cat);
}


void freeSourceCatalog(SourceCatalog** cat)
{
  if (NULL!=*cat) {
    // Free the KD-Tree.
    if (NULL!=(*cat)->tree) {
      freeKDTreeElement(&((*cat)->tree));
    }
    // Free the SIMPUT source catalog.
    if (NULL!=(*cat)->simput) {
      freeSimputSourceCatalog(&((*cat)->simput));
    }
    free(*cat);
    *cat=NULL;
  }    
}


SourceCatalog* loadSourceCatalog(const char* const filename,
				     const GenDet* const det,
				     int* const status)
{
  headas_chat(3, "load source catalog from file '%s' ...\n", filename);

  SourceCatalog* cat = newSourceCatalog(status);
  CHECK_STATUS_RET(*status, cat);

  // Set reference to ARF for SIMPUT library.
  simputSetARF(det->arf);

  // Set refernce to the random number generator to be used by the
  // SIMPUT library routines.
  simputSetRndGen(sixt_get_random_number);

  // Use the routines from the SIMPUT library to load the catalog.
  cat->simput = loadSimputSourceCatalog(filename, status);
  CHECK_STATUS_RET(*status, cat);

  // Allocate memory for an array of the sources, which will be
  // converted into a KDTree later.
  Source* list = (Source*)malloc(cat->simput->nentries*sizeof(Source));
  CHECK_NULL_RET(list, *status,
		 "memory allocation for source list failed", cat);

  // Empty template object.
  Source* templatesrc = newSource(status);
  CHECK_STATUS_RET(*status, cat);

  // Loop over all entries in the SIMPUT source catalog.
  int ii;
  for (ii=0; ii<cat->simput->nentries; ii++) {

    // Start with an empty Source object for each entry.
    list[ii] = *templatesrc;

    // Set the properties from the SIMPUT catalog.
    list[ii].src = cat->simput->entries[ii];
  } 
  CHECK_STATUS_RET(*status, cat);
  // END of loop over all entries in the FITS table.

  // Build a KDTree from the source list (array of Source objects).
  cat->tree = buildKDTree2(list, cat->simput->nentries, 0, status);
  CHECK_STATUS_RET(*status, cat);

  // In a later development stage this could be directly stored in 
  // a memory-mapped space.

  // Release memory.
  if (templatesrc) free(templatesrc);
  if (list) free(list);

  return(cat);
}


LinkedPhoListElement* genFoVXRayPhotons(SourceCatalog* const cat, 
					const Vector* const pointing, 
					const float fov,
					const double t0, const double t1,
					int* const status)
{
  // Minimum cos-value for point sources close to the FOV (in the direct
  // neighborhood).
  const double close_mult = 1.2; 
  const double close_fov_min_align = cos(close_mult*fov/2.); 

  // Perform a range search over all sources in the KDTree and 
  // generate new photons for the sources within the FoV.
  LinkedPhoListElement* list = 
    KDTreeRangeSearch(cat->tree, 0, pointing, close_fov_min_align, 
		      t0, t1, status);

  return(list);
}


