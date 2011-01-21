#include "xraysourcecatalog.h"


XRaySourceCatalog* newXRaySourceCatalog(int* const status)
{
  XRaySourceCatalog* cat = (XRaySourceCatalog*)malloc(sizeof(XRaySourceCatalog));
  CHECK_NULL(cat, *status,
	     "memory allocation for XRaySourceCatalog failed");

  // Initialize pointers with NULL.
  cat->tree = NULL;
  cat->spectra = NULL;

  // Initialize values.
  cat->nspectra = 0;

  return(cat);
}


void freeXRaySourceCatalog(XRaySourceCatalog** cat)
{
  if (NULL!=*cat) {
    if (NULL!=(*cat)->tree) {
      // TODO Free the binary tree.
    }

    // Free the spectra.
    if (NULL!=(*cat)->spectra) {
      long ii;
      for (ii=0; ii<(*cat)->nspectra; ii++) {
	freeXRaySourceSpectrum(&((*cat)->spectra[ii]));
      }
      free((*cat)->spectra);
    }
    free(*cat);
    *cat=NULL;
  }    
}


XRaySourceCatalog* loadSourceCatalog(const char* const filename,
				     const GenDet* const det,
				     int* const status)
{
  headas_chat(3, "load source catalog from file '%s' ...\n", filename);

  XRaySourceCatalog* cat = newXRaySourceCatalog(status);
  CHECK_STATUS_RET(*status, cat);

  // Open the FITS file with the source catalog (SRC_CAT) extension.
  fitsfile* fptr;
  fits_open_table(&fptr, filename, READONLY, status);
  CHECK_STATUS_RET(*status, cat);
  
  int hdunum, hdutype;
  // After opening the FITS file, get the number of the current HDU.
  if (1==fits_get_hdu_num(fptr, &hdunum)) {
    // This is the primary array, so try to move to the first extension 
    // and see if it is a table.
    fits_movabs_hdu(fptr, 2, &hdutype, status);
  } else {
    // Get the HDU type.
    fits_get_hdu_type(fptr, &hdutype, status);
  }
  CHECK_STATUS_RET(*status, cat);
  
  // If the current HDU is an image extension, throw an error message:
  if (IMAGE_HDU==hdutype) {
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "Error: FITS extension in file '%s' is not a table "
	    "but an image (HDU number: %d)\n", filename, hdunum);
    SIXT_ERROR(msg);
    return(cat);
  }

  // Determine the number of sources.
  long nrows;
  fits_get_num_rows(fptr, &nrows, status);
  CHECK_STATUS_RET(*status, cat);

  // Determine the column number in the FITS table.
  int cra, cdec, cemin, cemax, cflux, cspectrum;
  fits_get_colnum(fptr, CASEINSEN, "RA", &cra, status);
  fits_get_colnum(fptr, CASEINSEN, "DEC", &cdec, status);
  fits_get_colnum(fptr, CASEINSEN, "E_MIN", &cemin, status);
  fits_get_colnum(fptr, CASEINSEN, "E_MAX", &cemax, status);
  fits_get_colnum(fptr, CASEINSEN, "FLUX", &cflux, status);
  fits_get_colnum(fptr, CASEINSEN, "SPECTRUM", &cspectrum, status);
  CHECK_STATUS_RET(*status, cat);

  // Loop over all entries in the FITS table.
  long row;
  for (row=0; row<nrows; row++) {

    // New XRaySource object for each entry.
    XRaySource* src = newXRaySource(status);
    CHECK_STATUS_BREAK(*status);

    // Read the data from the FITS table.
    int anynul;
    fits_read_col(fptr, TFLOAT, cra, row+1, 1, 1, &src->ra, 
		  &src->ra, &anynul, status);
    fits_read_col(fptr, TFLOAT, cdec, row+1, 1, 1, &src->dec, 
		  &src->dec, &anynul, status);
    float flux=0., emin=0., emax=0.;
    fits_read_col(fptr, TFLOAT, cflux, row+1, 1, 1, &flux, 
		  &flux, &anynul, status);
    fits_read_col(fptr, TFLOAT, cemin, row+1, 1, 1, &emin, 
		  &emin, &anynul, status);
    fits_read_col(fptr, TFLOAT, cemax, row+1, 1, 1, &emax, 
		  &emax, &anynul, status);
    char spec_filename[MAXFILENAME];
    strcpy(spec_filename, "");
    fits_read_col(fptr, TSTRING, cspectrum, row+1, 1, 1, &spec_filename, 
		  &spec_filename, &anynul, status);
    CHECK_STATUS_BREAK(*status);
    
    // Load the required spectrum, if not available yet.
    long ii;
    for (ii=0; ii<cat->nspectra; ii++) {
      if (0==strcmp(cat->spectra[ii]->filename, spec_filename)) break;
    }
    if (ii==cat->nspectra) {
      // The spectrum does not exist in the catalog yet.
      cat->spectra = realloc(cat->spectra, 
			     (cat->nspectra+1)*sizeof(XRaySourceSpectrum*));
      CHECK_NULL_BREAK(cat->spectra, *status,
		       "memory allocation for spectrum catalog failed");
      cat->nspectra++;
      cat->spectra[ii] = loadXRaySpectrumFilename(spec_filename, status);
      CHECK_STATUS_BREAK(*status);

      // Apply the ARF.
      applyARF2Spectrum(cat->spectra[ii], det->arf, status);
      CHECK_STATUS_BREAK(*status);      
    }
    // Assign the spectrum to the source.
    src->nspectra   = 1;
    src->spectra    = (XRaySourceSpectrum**)malloc(sizeof(XRaySourceSpectrum*));
    CHECK_NULL_BREAK(cat->spectra, *status,
		     "memory allocation for spectrum list failed");    
    src->spectra[0] = cat->spectra[ii];

    // Determine the photon rate from this particular source.
    src->pps = 
      getSpectralPhotonRate(src->spectra[0], emin, emax) *
      flux/getSpectralEnergyFlux(src->spectra[0], emin, emax);

    // Insert the XRaySource object in the KDTree.
    // TODO
    
  } 
  CHECK_STATUS_RET(*status, cat);
  // END of loop over all entries in the FITS table.

  // Close the FITS file.
  fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, cat);

  return(cat);
}

