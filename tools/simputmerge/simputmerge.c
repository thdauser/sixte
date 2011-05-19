#include "simputmerge.h"


int simputmerge_main() 
{
  // Program parameters.
  struct Parameters par;

  // Filenames of the input catalogs.
  char infilenames[2][MAXFILENAME];

  // SIMPUT source catalogs.
  SimputSourceCatalog* incat[2]={NULL, NULL};
  SimputSourceCatalog* outcat  = NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL
  set_toolname("simputmerge");
  set_toolversion("0.01");


  do { // Beginning of ERROR HANDLING Loop.

    // ---- Initialization ----
    
    // Read the parameters using PIL.
    status=simputmerge_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Load the input catalogs.
    strcpy(infilenames[0], par.Infile1);
    strcpy(infilenames[1], par.Infile2);
    int ii;
    for (ii=0; ii<2; ii++) {
      incat[ii] = loadSimputSourceCatalog(infilenames[ii], &status);
      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);

    // Get an empty object for the output catalog.
    outcat = getSimputSourceCatalog(&status);
    CHECK_STATUS_BREAK(status);

    // Allocate sufficient memory to store the entries of both
    // input catalogs in the output catalog.
    outcat->entries = 
      (SimputSourceEntry**)
      malloc((incat[0]->nentries+incat[1]->nentries)*sizeof(SimputSourceEntry*));
    CHECK_NULL_BREAK(outcat->entries, status, "memory allocation failed!\n");

    // Smallest still available SRC_ID.
    long min_src_id=1; 
    // Loop over both source catalogs.
    for (ii=0; ii<2; ii++) {
      // Loop over all entries in the source catalog.
      long jj;
      for (jj=0; jj<incat[ii]->nentries; jj++) {

	// Check if the SRC_ID of the new source is already contained 
	// in the output catalog. The SRC_ID entry must be unique.
	long src_id = incat[ii]->entries[jj]->src_id;
	if (src_id<min_src_id) {
	  src_id=min_src_id;
	}
	long kk;
	do {
	  for (kk=0; kk<outcat->nentries; kk++) {
	    if (src_id==outcat->entries[kk]->src_id) {
	      // This SRC_ID is not available any more.
	      if (src_id==min_src_id) {
		min_src_id++;
	      }
	      src_id++;
	      break;
	    }
	  } 
	} while (kk<outcat->nentries);
	if (src_id==min_src_id) {
	  min_src_id++;
	}

	// Handle spectrum, image, and lightcur extensions.
	char spectrum[MAXFILENAME];
	char lightcur[MAXFILENAME];
	char image[MAXFILENAME];

	// Check whether the extensions should remain in their current
	// place or if they should by copied to the new output file.
	if (0==par.FetchExtensions) {
	  // Extensions should remain in their original place.
	  // Check if they are local references to the same file
	  // containing the catalog.
	  // Spectrum.
	  if ('['==incat[ii]->entries[jj]->spectrum[0]) {
	    strcpy(spectrum, infilenames[ii]);
	    strcat(spectrum, incat[ii]->entries[jj]->spectrum);
	  } else {
	    strcpy(spectrum, incat[ii]->entries[jj]->spectrum);
	  }
	  // Image.
	  if ('['==incat[ii]->entries[jj]->image[0]) {
	    strcpy(image, infilenames[ii]);
	    strcat(image, incat[ii]->entries[jj]->image);
	  } else {
	    strcpy(image, incat[ii]->entries[jj]->image);
	  }
	  // Light curve.
	  if ('['==incat[ii]->entries[jj]->lightcur[0]) {
	    strcpy(lightcur, infilenames[ii]);
	    strcat(lightcur, incat[ii]->entries[jj]->lightcur);
	  } else {
	    strcpy(lightcur, incat[ii]->entries[jj]->lightcur);
	  }
	  // TODO What happens if the input file resides in another 
	  // directory than the output file.
	} else {
	  // Extensions should be copied to the new output file.

	  // TODO
	}	
	// END of extensions should be copied to the new output file.

	// Copy the entry from the input to the output catalog.
	outcat->nentries++;
	outcat->entries[outcat->nentries-1] = 
	  getSimputSourceEntryV(src_id, 
				incat[ii]->entries[jj]->src_name,
				incat[ii]->entries[jj]->ra,
				incat[ii]->entries[jj]->dec,
				incat[ii]->entries[jj]->imgrota,
				incat[ii]->entries[jj]->imgscal,
				incat[ii]->entries[jj]->e_min,
				incat[ii]->entries[jj]->e_max,
				incat[ii]->entries[jj]->flux,
				spectrum, image, lightcur,
				&status);
	CHECK_STATUS_BREAK(status);

	// Output of progress.
	if (0==outcat->nentries % 100) {
	  headas_chat(1, "\r%ld/%ld (%.1lf%%) entries", 
		      outcat->nentries, incat[0]->nentries+incat[1]->nentries,
		      outcat->nentries*100./(incat[0]->nentries+incat[1]->nentries));
	  fflush(NULL);
	}
      }
      CHECK_STATUS_BREAK(status);
      // END of loop over all entries in the source catalog.
    }
    CHECK_STATUS_BREAK(status);
    headas_chat(1, "\n");
    // END of loop over both source catalogs.

    // Store the output catalog in the FITS file.
    saveSimputSourceCatalog(outcat, par.Outfile, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of ERROR HANDLING Loop.

  // --- Clean up ---
  headas_chat(3, "\ncleaning up ...\n");

  // Release memory.
  freeSimputSourceCatalog(&incat[0]);
  freeSimputSourceCatalog(&incat[1]);
  freeSimputSourceCatalog(&outcat);

  if (status==EXIT_SUCCESS) headas_chat(0, "finished successfully!\n\n");
  return(status);
}



int simputmerge_getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status = EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("Infile1", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the input file 1!\n", status);
    return(status);
  } 
  strcpy(par->Infile1, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("Infile2", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the input file 2!\n", status);
    return(status);
  } 
  strcpy(par->Infile2, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("Outfile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the output file!\n", status);
    return(status);
  } 
  strcpy(par->Outfile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_bool("FetchExtensions", &par->FetchExtensions);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the FetchExtensions parameter!\n", status);
    return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the clobber parameter!\n", status);
    return(status);
  }

  return(status);
}


