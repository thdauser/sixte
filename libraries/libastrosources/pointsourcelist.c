#include "pointsourcelist.h"


PointSourceList* selectFoVPointSourceList(kdNode* tree, 
					  PointSourceFile* psf,
					  Vector* telescope_direction,
					  const double max_align,
					  int* status)
{
  // Allocate memory for a new PointSourceList data structure.
  PointSourceList* psl = (PointSourceList*)malloc(sizeof(PointSourceList));
  if (NULL==psl) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: memory allocation for PointSourceList failed!\n",
		   *status);
    return(psl);
  }
  // Set default values:
  psl->nsources = 0;
  psl->sources  = NULL;
  
  // Perform a range search on the kd-Tree and add the selected
  // sources to a new PointSourceList combining the source positions
  // with the other important data like the source count rate or 
  // spectrum.
  SourceList* sl=NULL;
  long nelements=0;
  *status=kdTreeRangeSearch(tree, 0, telescope_direction, 
			    max_align*max_align, 
			    &sl, &nelements);
  if (EXIT_SUCCESS!=*status) return(psl);

  // Allocate memory for the entries PointSourceList.
  psl->nsources = nelements;
  psl->sources  = (PointSource*)malloc(nelements*sizeof(PointSource));

  // Add the sources in the SourceList to the PointSourceList 
  // and associate the data from the PointSourceFile.
  PointSource ps; // Buffer for reading from the PointSourceFile.
  long count;
  for(count=0; count<nelements; count++) {

    // Read the source data from the PointSourceFile.
    get_PointSourceTable_Row(psf, sl[count].line, &ps, status);

    // Set the spectrum and light curve of the PointSource:
    // Source spectrum.
    if ((ps.spectrum_index<1) || 
	(ps.spectrum_index>psf->spectrumstore.nspectra)) {
      headas_chat(0, "\n### Warning: no source spectrum specified for point source!\n"
		  "     Using default spectrum instead.\n");
      ps.spectrum = &psf->spectrumstore.spectrum[0];
    } else {
      ps.spectrum = &psf->spectrumstore.spectrum[ps.spectrum_index-1];
    } 
    // END of assigning a spectrum to the source.

    // Light curve.
    if (0<ps.lc_type) {
      // Try and find out the name of the FITS file containing the light curve.
      char lc_filename[MAXMSG], key[MAXMSG], comment[MAXMSG]; // Buffers
      sprintf(key, "LC%06ld", ps.lc_type);
      if (fits_read_key(psf->fptr, TSTRING, key, 
			lc_filename, comment, status)) break;
      // Then load the file data and store the light curve.
      ps.lc = loadLinLightCurveFromFile(lc_filename, ps.rate, status);
      if (EXIT_SUCCESS!=*status) break;
      // So far there was no photon created for this source.
      ps.t_last_photon = 0.;
	  
    } else {
      // No particular light curve has been assigned to this source.
      // Set light curve pointer to NULL.
      ps.lc = NULL;
      // So far there was no photon created for this source.
      ps.t_last_photon = 0.;
	  
    } // END of assigning a light curve to the source.

    // Assign the data to the entry in the PointSourceList.
    psl->sources[count] = ps;

  }
  // END of loop over all sources in the SourceList.

  freeSourceList(sl);

  return(psl);
}



void clearPointSourceList(PointSourceList* psl)
{
  if (NULL!=psl) {
    if (NULL!=psl->sources) {

      // Free the light curves of the individual sources.
      long count;
      for (count=0; count<psl->nsources; count++) {
	if (psl->sources[count].lc != NULL) {
	  freeLinLightCurve(psl->sources[count].lc);
	  psl->sources[count].lc=NULL;
	}
      }

      free(psl->sources);
      psl->nsources = 0;
    }
    free(psl);
  }
}


