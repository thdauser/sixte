#include "sourcecatalog.h"


SourceCatalog* newSourceCatalog(int* const status)
{
  SourceCatalog* cat = (SourceCatalog*)malloc(sizeof(SourceCatalog));
  CHECK_NULL(cat, *status,
	     "memory allocation for SourceCatalog failed");

  // Initialize pointers with NULL.
  cat->tree       =NULL;
  cat->extsources =NULL;
  cat->nextsources=0;
  cat->simput     =NULL;

  return(cat);
}


void freeSourceCatalog(SourceCatalog** cat)
{
  if (NULL!=*cat) {
    // Free the KD-Tree.
    if (NULL!=(*cat)->tree) {
      freeKDTreeElement(&((*cat)->tree));
    }
    // Free the array of extended sources.
    if (NULL!=(*cat)->extsources) {
      free((*cat)->extsources);
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

  // Determine the number of point-like and the number of 
  // extended sources.
  long nextended  = 0;
  long npointlike = 0;
  long ii;
  for (ii=0; ii<cat->simput->nentries; ii++) {
    float extension = getSimputSourceExtension(cat->simput->entries[ii], status);
    CHECK_STATUS_BREAK(*status);
    if (extension>0.) {
      // This is an extended source.
      nextended++;
    } else {
      // This is a point-like source.
      npointlike++;
    }
  }
  CHECK_STATUS_RET(*status, cat);

  // Allocate memory for an array of the point-like sources, 
  // which will be converted into a KDTree afterwards.
  Source* list = (Source*)malloc(npointlike*sizeof(Source));
  CHECK_NULL_RET(list, *status,
		 "memory allocation for source list failed", cat);

  // Allocate memory for the array of the extended sources.
  cat->extsources = (Source*)malloc(cat->nextsources*sizeof(Source));
  CHECK_NULL_RET(list, *status,
		 "memory allocation for source list failed", cat);

  // Empty template object.
  Source* templatesrc = newSource(status);
  CHECK_STATUS_RET(*status, cat);

  // Loop over all entries in the SIMPUT source catalog.
  long cpointlike =0;
  cat->nextsources=0;
  for (ii=0; ii<cat->simput->nentries; ii++) {
    float extension = getSimputSourceExtension(cat->simput->entries[ii], status);
    CHECK_STATUS_BREAK(*status);
    if (extension>0.) {
      // This is an extended source.
      cat->nextsources++;
      // Start with an empty Source object for this entry.
      cat->extsources[cat->nextsources-1] = *templatesrc;
      // Set the properties from the SIMPUT catalog.
      cat->extsources[cat->nextsources-1].src = cat->simput->entries[ii];

    } else {
      // This is a point-like source.
      cpointlike++;
      list[cpointlike-1]     = *templatesrc;
      list[cpointlike-1].src = cat->simput->entries[ii];
    }
  } 
  CHECK_STATUS_RET(*status, cat);
  // END of loop over all entries in the FITS table.

  // Build a KDTree from the source list (array of Source objects).
  cat->tree = buildKDTree2(list, npointlike, 0, status);
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
  // The kdTree only contains point-like sources.
  LinkedPhoListElement* list = 
    KDTreeRangeSearch(cat->tree, 0, pointing, close_fov_min_align, 
		      t0, t1, status);

  // Loop over all extended sources.
  long ii;
  for (ii=0; ii<cat->nextsources; ii++) {
    // Get the maximum extension of the source.
    float extension=
      getSimputSourceExtension(cat->extsources[ii].src, status);
    CHECK_STATUS_RET(*status, list);

    // Check if at least a part of the source lies within the FoV.
    Vector location=unit_vector(cat->extsources[ii].src->ra, 
				cat->extsources[ii].src->dec);
    if (fabs(scalar_product(&location, pointing)) > 
	cos(close_mult*(fov*0.5 + extension))) {
      // Generate photons for this particular source.
      LinkedPhoListElement* newlist = 
	getXRayPhotons(&(cat->extsources[ii]), t0, t1, status);
      CHECK_STATUS_RET(*status, list);

      // Merge the new photons into the existing list.
      list = mergeLinkedPhoLists(list, newlist);
    }
  }

  return(list);
}


