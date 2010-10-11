#include "pointsourcecatalog.h"


PointSourceCatalog openPointSourceCatalog(char* filename, int hdu, int* status)
{
  PointSourceCatalog psc;
  
  // Set default initial values.
#ifdef POINTSOURCE_KDTREE
  psc.kdtree = NULL;
#else
  psc.psl = NULL;
#endif

  // Open the PointSourceFile.
  psc.file = openPointSourceFile(filename, hdu, status);

  return(psc);
}



void clearPointSourceCatalog(PointSourceCatalog* psc)
{
  if (NULL!=psc) {
#ifdef POINTSOURCE_KDTREE
    freeKDTree(psc->kdtree);
    psc->kdtree=NULL;
#else
    freePointSourceList(psc->psl);
    psc->psl=NULL;
#endif
    free_PointSourceFile(psc->file);
  }
}



void preselectPointSources(PointSourceCatalog* psc, Vector normal, 
			   double max_align, int* status)
{
  // Allocate memory for the PointSourceList.
  PointSourceList* psl=NULL;
  psl = (PointSourceList*)malloc(sizeof(PointSourceList));
  if (NULL==psl) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for PointSourceList failed!\n",
		   *status);
    return;
  }
  psl->nsources = 0;
  psl->sources = (PointSource*)malloc(MAX_N_POINTSOURCES*sizeof(PointSource));
  if (NULL==psl->sources) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for PointSourceList failed!\n",
		   *status);
    return;
  }

  // Read all sources from the FITS file.
  long row;
  PointSource* ps=psl->sources; // Pointer to the last entry in PointSourceList.

  for (row=0; row<psc->file->nrows; row++) {
    
    // Read the PointSource from the FITS table.
    get_PointSourceTable_Row(psc->file, row, ps, status);
    if (EXIT_SUCCESS!=*status) return;

    // Check, whether the PointSource is within the preselection 
    // band.
    if(fabs(scalar_product(&ps->location, &normal)) < max_align) {
      if(psl->nsources+1 >= MAX_N_POINTSOURCES) {
	// Too many sources in the PointSourceCatalog!
	*status=EXIT_FAILURE;
	char msg[MAXMSG]; // Error output buffer.
	sprintf(msg, "Error: too many sources (%ld) in the PointSourceCatalog!\n", 
		psl->nsources+1);
	HD_ERROR_THROW(msg, *status);
	return;
      }

      // If the light curve type is a value > 0 meaning that the particular
      // light curve for this source is given in a FITS file, load that file
      // and assign the pointer to the light curve.
      if (0<ps->lc_type) {
	// Try and find out the name of the FITS file containing the light curve.
	char lc_filename[MAXMSG], key[MAXMSG], comment[MAXMSG]; // Buffers
	sprintf(key, "LC%06ld", ps->lc_type);
	if (fits_read_key(psc->file->fptr, TSTRING, key, 
			  lc_filename, comment, status)) return;
	// Then load the file data and store the light curve.
	ps->lc = loadLinLightCurveFromFile(lc_filename, ps->rate, 
					   status);
	if (EXIT_SUCCESS!=*status) return;
	// So far there was no photon created for this source.
	ps->t_last_photon = 0.;
	
      } else {
	// No particular light curve has been assigned to this source.
	// Set light curve pointer to NULL.
	ps->lc = NULL;
	// So far there was no photon created for this source.
	ps->t_last_photon = 0.;
	  
      } 
      // END of assigning a light curve to the source.
	  
      // Source spectrum.
      if ((ps->spectrum_index<1) || 
	  (ps->spectrum_index>psc->file->spectrumstore.nspectra)) {
	headas_chat(0, "\n### Warning: no source spectrum specified for point source!\n"
		    "     Using default spectrum instead.\n");
	ps->spectrum = &psc->file->spectrumstore.spectrum[0];
      } else {
	ps->spectrum = &psc->file->spectrumstore.spectrum[ps->spectrum_index-1];
      } 
      // END of assigning a spectrum to the source.
      
      // Increase number of sources in the selected catalog
      psl->nsources++;
      ps = &psl->sources[psl->nsources];
    } 
    // END of check whether this sources lies within the preselection band.
  }
  // END of reading all sources from the FITS file.

#ifdef POINTSOURCE_KDTREE
  // TODO RM
  printf("\npreselection with %ld sources ...\n", psl->nsources);

  // Destroy the old KDTree.
  freeKDTree(psc->kdtree);

  // Build a kdTree from the preselected PointSourceList.
  if (psl->nsources>0) {
    psc->kdtree = buildKDTree(psl->sources, psl->nsources);
    if (NULL==psc->kdtree) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: Building kd-Tree failed!\n", *status);
      return;
    }
  } else {
    psc->kdtree=NULL;
  }

  // Release the memory of the PointSourceList.
  free(psl->sources);
  psl->sources=NULL;
  psl->nsources=0;
  free(psl);

#else 

  freePointSourceList(psc->psl);
  psc->psl = psl;

#endif
}



void generateFoVPointSourcePhotons(PointSourceCatalog* psc,
				   Vector* ref, double min_align, 
				   double time, double dt, 
				   struct PhotonOrderedListEntry** pl,
				   const struct ARF* const arf,
				   int* status)
{
#ifdef POINTSOURCE_KDTREE
  kdnfound=0; // RM
  kdnchecked=0;
  kdTreeRangeSearch(psc->kdtree, 0, ref, min_align, 
		    time, dt, pl, arf, status);
#else

  long count;
  for (count=0; count<psc->psl->nsources; count++) {
    // Check if the source is close to / within the FoV.
    if(fabs(scalar_product(&psc->psl->sources[count].location, ref)) 
       > min_align) {

      *status=create_PointSourcePhotons(&psc->psl->sources[count],
					time, dt, pl, arf);
      if (EXIT_SUCCESS!=*status) return;
    }
    // END of check if source is close to the FoV.
  }
  // END of loop over all sources in the preselected 
  // PointSourceList.

#endif
}



