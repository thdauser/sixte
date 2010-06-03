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
  PointSource ps; // Input buffer.
  // Unit vector pointing in the direction of the source.
  Vector source_direction;

  for (row=0; row<psc->file->nrows; row++) {

    // Read the PointSource from the FITS table.
    get_PointSourceTable_Row(psc->file, row, &ps, status);
    if (EXIT_SUCCESS!=*status) return;

    // Get a unit vector pointing in the direction of the source:
    source_direction = unit_vector(ps.ra, ps.dec);

    // Check, whether the PointSource is within the preselection 
    // band.
    if(fabs(scalar_product(&source_direction, &normal)) < max_align) {
      if(psl->nsources+1 > MAX_N_POINTSOURCES) {
	// Too many sources in the PointSourceCatalog !
	*status=EXIT_FAILURE;
	char msg[MAXMSG];  // Error output buffer.
	sprintf(msg, "Error: too many sources (%ld) in the PointSourceCatalog!\n", 
		psl->nsources+1);
	HD_ERROR_THROW(msg, *status);
	return;
      }

      // Add the current source to the selected catalog:
      psl->sources[psl->nsources].rate = ps.rate;
      
      // Save the source direction in the source catalog-array:
      psl->sources[psl->nsources].ra  = ps.ra;
      psl->sources[psl->nsources].dec = ps.dec;

      // Set the light curve type to the value specified in the source catalog.
      psl->sources[psl->nsources].lc_type = ps.lc_type;
      // If the light curve type is a value > 0 meaning that the particular
      // light curve for this source is given in a FITS file, load that file
      // and assign the pointer to the light curve.
      if (0<ps.lc_type) {
	// Try and find out the name of the FITS file containing the light curve.
	char lc_filename[MAXMSG], key[MAXMSG], comment[MAXMSG]; // Buffers
	sprintf(key, "LC%06ld", ps.lc_type);
	if (fits_read_key(psc->file->fptr, TSTRING, key, 
			  lc_filename, comment, status)) return;
	// Then load the file data and store the light curve.
	psl->sources[psl->nsources].lc = 
	  loadLinLightCurveFromFile(lc_filename, ps.rate, status);
	if (EXIT_SUCCESS!=*status) return;
	// So far there was no photon created for this source.
	psl->sources[psl->nsources].t_last_photon = 0.;
	
      } else {
	// No particular light curve has been assigned to this source.
	// Set light curve pointer to NULL.
	psl->sources[psl->nsources].lc = NULL;
	// So far there was no photon created for this source.
	psl->sources[psl->nsources].t_last_photon = 0.;
	  
      } 
      // END of assigning a light curve to the source.
	  
      // Source spectrum.
      if ((ps.spectrum_index<1) || 
	  (ps.spectrum_index>psc->file->spectrumstore.nspectra)) {
	headas_chat(0, "\n### Warning: no source spectrum specified for point source!\n"
		    "     Using default spectrum instead.\n");
	psl->sources[psl->nsources].spectrum = &psc->file->spectrumstore.spectrum[0];
      } else {
	psl->sources[psl->nsources].spectrum = 
	  &psc->file->spectrumstore.spectrum[ps.spectrum_index-1];
      } 
      // END of assigning a spectrum to the source.
      
      // Increase number of sources in the selected catalog
      psl->nsources++;
    } 
    // END of check whether this sources lies within the preselection band.
  }
  // END of reading all sources from the FITS file.

#ifdef POINTSOURCE_KDTREE
  // Build a kdTree from the preselected PointSourceList.
  psc->kdtree = buildKDTree(psl->sources, psl->nsources);
  if (NULL==psc->kdtree) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Building kd-Tree failed!\n", *status);
    return;
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



LinkedPointSourceList getFoVPointSources(PointSourceCatalog* psc,
					 Vector* ref, double min_align, 
					 int* status)
{
#ifdef POINTSOURCE_KDTREE
  return(kdTreeRangeSearch(psc->kdtree, 0, ref, pow(min_align,2.), status));
#else

  LinkedPointSourceList lpsl=NULL;
  LinkedPointSourceListEntry** entry=&lpsl;

  // Loop over all sources in the preselected catalog.
  long count;
  Vector source_direction;
  for(count=0; count<psc->psl->nsources; count++) {
    // Check if the source is close to the FoV.
    source_direction = unit_vector(psc->psl->sources[count].ra,
				   psc->psl->sources[count].dec);
    if(fabs(scalar_product(&source_direction, ref)) > min_align) {
      // Add the source to the LinkedPointSourceList.
      *entry = (LinkedPointSourceListEntry*)malloc(sizeof(LinkedPointSourceListEntry));
      if (NULL==*entry) {
	*status=EXIT_FAILURE;
	HD_ERROR_THROW("Error: Memory allocation for LinkedPointSourceList failed!\n",
		       *status);
	return(lpsl);
      }
      // Check if this is the first entry in the list.
      if (NULL==lpsl) {
	lpsl = (*entry);
      }
      (*entry)->source = &psc->psl->sources[count];
      entry = &(*entry)->next;
    }
    // End of check whether the source is within / close to the FoV.
  }

  return(lpsl);
#endif
}



// OBSOLETE
/*
PointSourceCatalog* getPointSourceCatalog(PointSourceFile* psf, 
					  Vector normal_vector, 
					  const double max_align,
					  int* status)
{
  PointSourceCatalog* psc = NULL;
  char msg[MAXMSG];  // Error output buffer.
  
  do { // Beginning of outer ERROR handling loop.
    
    // Allocate memory:
    psc = (PointSourceCatalog*)malloc(sizeof(PointSourceCatalog));
    if (NULL==psc) {
      *status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: not enough memory available to store the PointSourceCatalog!\n", 
		     *status);
      break;
    }
    psc->nsources = 0;
    
    psc->sources = (PointSource*)malloc(MAX_N_POINTSOURCES*sizeof(PointSource));
    if(NULL==psc->sources) {
      *status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: not enough memory available to store the PointSourceCcatalog!\n", 
		     *status);
      break;
    } 
    // END of memory allocation


    // At the moment the selected catalog is empty.
    PointSource ps;
    // Counter to access the individual sources in the catalog.
    int source_counter; 
    // Unit vector pointing in the direction of the source.
    Vector source_direction;

    // Read-in all sources from the individual table rows, one after another:
    for (source_counter=0; source_counter < psf->nrows; source_counter++) {
      
      // Read source data (right asension, declination, photon rate, spectrum):
      if (get_PointSourceTable_Row(psf, source_counter, &ps, status)) break;

      // Get a unit vector pointing in the direction of the source:
      source_direction = unit_vector(ps.ra, ps.dec);
	  
      // Check whether the source should be added to the preselected catalog:
      if(fabs(scalar_product(&source_direction, &normal_vector)) < max_align) {
	if(psc->nsources > MAX_N_POINTSOURCES) {
	  // Too many sources in the PointSourceCatalog !
	  *status=EXIT_FAILURE;
	  sprintf(msg, "Error: too many sources (%ld) in the PointSourceCatalog!\n", 
		  psc->nsources+1);
	  HD_ERROR_THROW(msg, *status);
	  break;
	}

	// Add the current source to the selected catalog:
	psc->sources[psc->nsources].rate = ps.rate;

	// Save the source direction in the source catalog-array:
	psc->sources[psc->nsources].ra  = ps.ra;
	psc->sources[psc->nsources].dec = ps.dec;

	// Set the light curve type to the value specified in the source catalog.
	psc->sources[psc->nsources].lc_type = ps.lc_type;
	// If the light curve type is a value > 0 meaning that the particular
	// light curve for this source is given in a FITS file, load that file
	// and assign the pointer to the light curve.
	if (0<ps.lc_type) {
	  // Try and find out the name of the FITS file containing the light curve.
	  char lc_filename[MAXMSG], key[MAXMSG], comment[MAXMSG]; // Buffers
	  sprintf(key, "LC%06ld", ps.lc_type);
	  if (fits_read_key(psf->fptr, TSTRING, key, 
			    lc_filename, comment, status)) break;
	  // Then load the file data and store the light curve.
	  psc->sources[psc->nsources].lc = 
	    loadLinLightCurveFromFile(lc_filename, ps.rate, status);
	  if (EXIT_SUCCESS!=*status) break;
	  // So far there was no photon created for this source.
	  psc->sources[psc->nsources].t_last_photon = 0.;
	  
	} else {
	  // No particular light curve has been assigned to this source.
	  // Set light curve pointer to NULL.
	  psc->sources[psc->nsources].lc = NULL;
	  // So far there was no photon created for this source.
	  psc->sources[psc->nsources].t_last_photon = 0.;
	  
	} // END of assigning a light curve to the source.
	  
	// Source spectrum.
	if ((ps.spectrum_index<1) || 
	    (ps.spectrum_index>psf->spectrumstore.nspectra)) {
	  headas_chat(0, "\n### Warning: no source spectrum specified for point source!\n"
		      "     Using default spectrum instead.\n");
	  psc->sources[psc->nsources].spectrum = &psf->spectrumstore.spectrum[0];
	} else {
	  psc->sources[psc->nsources].spectrum = 
	    &psf->spectrumstore.spectrum[ps.spectrum_index-1];
	} // END of assigning a spectrum to the source.
	
	// Increase number of sources in the selected catalog
	psc->nsources++;
      } // END of check whether this sources lies within the preselection band.
    } // END of loop over all sources in the current catalog file
    if (EXIT_SUCCESS!=*status) break;
  } while (0); // END of outer error handling loop

  return(psc);
}



void free_PointSourceCatalog(PointSourceCatalog* psc)
{
  if (NULL!=psc) {
    if (NULL!=psc->sources) {
      int count;
      for (count=0; count<psc->nsources; count++) {
	if (psc->sources[count].lc != NULL) {
	  freeLinLightCurve(psc->sources[count].lc);
	  psc->sources[count].lc=NULL;
	}
      }
      free(psc->sources);
      psc->sources=NULL;
    }
    free(psc);
  }
}
*/
