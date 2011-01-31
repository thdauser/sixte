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
    // Free the KD-Tree.
    if (NULL!=(*cat)->tree) {
      freeKDTreeElement(&((*cat)->tree));
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

  // Allocate memory for an array of the sources, which will be
  // converted into a KDTree later.
  XRaySource* list = (XRaySource*)malloc(nrows*sizeof(XRaySource));
  CHECK_NULL_RET(list, *status,
		 "memory allocation for source list failed", cat);

  // Empty template object.
  XRaySource* templatesrc = newXRaySource(status);
  CHECK_STATUS_RET(*status, cat);

  // Loop over all entries in the FITS table.
  long row;
  for (row=0; row<nrows; row++) {

    // Start with an empty XRaySource object for each entry.
    list[row] = *templatesrc;

    // Read the data from the FITS table.
    int anynul=0;
    fits_read_col(fptr, TFLOAT, cra, row+1, 1, 1, &(list[row].ra), 
		  &(list[row].ra), &anynul, status);
    fits_read_col(fptr, TFLOAT, cdec, row+1, 1, 1, &(list[row].dec), 
		  &(list[row].dec), &anynul, status);
    float flux=0., emin=0., emax=0.;
    fits_read_col(fptr, TFLOAT, cflux, row+1, 1, 1, &flux, 
		  &flux, &anynul, status);
    fits_read_col(fptr, TFLOAT, cemin, row+1, 1, 1, &emin, 
		  &emin, &anynul, status);
    fits_read_col(fptr, TFLOAT, cemax, row+1, 1, 1, &emax, 
		  &emax, &anynul, status);
    char buffer[MAXFILENAME]="";
    char* spec_filename = &(buffer[0]);
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
    list[row].nspectra = 1;
    list[row].spectra  = (XRaySourceSpectrum**)malloc(sizeof(XRaySourceSpectrum*));
    CHECK_NULL_BREAK(cat->spectra, *status,
		     "memory allocation for spectrum list failed");    
    list[row].spectra[0] = cat->spectra[ii];

    // Determine the photon rate from this particular source.
    list[row].pps = 
      getSpectralPhotonRate(list[row].spectra[0], emin, emax) *
      flux/getSpectralEnergyFlux(list[row].spectra[0], emin, emax);

  } 
  CHECK_STATUS_RET(*status, cat);
  // END of loop over all entries in the FITS table.

  // Build a KDTree from the source list (array of XRaySource objects).
  cat->tree = buildKDTree2(list, nrows, 0, status);
  CHECK_STATUS_RET(*status, cat);

  // In a later development stage this could be directly stored in 
  // a memory-mapped space.

  // Release memory.
  if (templatesrc) free(templatesrc);
  if (list) free(list);

  // Close the FITS file.
  fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, cat);

  return(cat);
}



LinkedPhoListElement* genFoVXRayPhotons(XRaySourceCatalog* const cat, 
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


