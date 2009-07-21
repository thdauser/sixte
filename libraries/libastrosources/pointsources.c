#include "pointsources.h"


////////////////////////////////////////////////////
PointSourceFileCatalog* get_PointSourceFileCatalog() {
  PointSourceFileCatalog* psfc = NULL;

  psfc = (PointSourceFileCatalog*) malloc(sizeof(PointSourceFileCatalog));

  if(psfc!=NULL) {
    psfc->nfiles = 0;
    psfc->files = NULL;
  }

  return(psfc);
}




/////////////////////////////////////////////////////////
void free_PointSourceFileCatalog(PointSourceFileCatalog* psfc) {
  if (NULL!=psfc) {
    if (psfc->nfiles>0) {
      int count;
      for(count=0; count<psfc->nfiles; count++) {
	free(psfc->files[count]);
      }
    }

    if (NULL != psfc->files) {
      free(psfc->files);
    }

    free(psfc);
  }
}




//////////////////////////////////////////////////////////
PointSourceFile* get_PointSourceFile() {
  PointSourceFile* psf = (PointSourceFile*)malloc(sizeof(PointSourceFile));

  if (NULL!=psf) {
    psf->fptr = NULL;
    psf->cra = 0;
    psf->cdec = 0;
    psf->crate = 0;
    psf->cspectrum = 0;
  }

  return(psf);
}


//////////////////////////////////////////////////////////
PointSourceFile* get_PointSourceFile_fromFile(char* filename, int* status) 
{
  PointSourceFile* psf = NULL;
  char msg[MAXMSG];      // error output buffer

  
  do { // Beginning of ERROR handling loop.
    
    // Call the generic constructor.
    psf = get_PointSourceFile();
    if (NULL==psf) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: memory allocation for PointSourceFile failed!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }

    // OPEN the specified FITS file and store basic information.
    headas_chat(5, "open point-source catalog in file '%s' ...\n", filename);

    // Open the source catalogue (FITS-file):
    if(fits_open_table(&psf->fptr, filename, READONLY, status)) break;
    // Get the HDU number in the source catalogue:
    int hdunum, hdutype; // Needed to access the HDUs in the source FITS file.
    if(1 == fits_get_hdu_num(psf->fptr, &hdunum)) {
      // This is the primary array,
      // so try to move to the first extension and see if it is a table:
      fits_movabs_hdu(psf->fptr, 2, &hdutype, status);
    } else {
      // Get the HDU type:
      fits_get_hdu_type(psf->fptr, &hdutype, status);
    }
    // Image HDU results in an error message:
    if (IMAGE_HDU==hdutype) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: FITS extension in source catalog file '%s' is "
	      "not a table but an image (HDU number: %d)!\n", filename, hdunum);
      HD_ERROR_THROW(msg, *status);
      break;
    }

    // Determine the column numbers of the right ascension, declination,
    // photon rate, and spectrum columns in the FITS table
    if (fits_get_colnum(psf->fptr, CASEINSEN, "RA", &psf->cra, status)) break;
    if (fits_get_colnum(psf->fptr, CASEINSEN, "DEC", &psf->cdec, status)) break;
    if (fits_get_colnum(psf->fptr, CASEINSEN, "PPS", &psf->crate, status)) break;
    if (fits_get_colnum(psf->fptr, CASEINSEN, "SPECTRUM", &psf->cspectrum, status)) break;

    // TODO: If spectrum column is not available, use HEADER keyword instead.

    // Determine the number of rows in the FITS table:
    if (fits_get_num_rows(psf->fptr, &psf->nrows, status)) break;	

  } while (0); // END of ERROR handling loop.

  return(psf);
}



///////////////////////////////////////////////////////////
void free_PointSourceFile(PointSourceFile* psf) {
  if (NULL!=psf) {

    // Close the FITS file if still open.
    if (NULL!=psf->fptr) {
      int status;
      fits_close_file(psf->fptr, &status);
    }

    free(psf);
  }
}




//////////////////////////////////////////////////////////////////////////
PointSourceCatalog* get_PointSourceCatalog(PointSourceFileCatalog* psfc, 
					   Vector normal_vector, 
					   const double max_align,
					   struct Spectrum_Store spectrum_store,
					   int* status)
{
  PointSourceCatalog* psc = NULL;

  char msg[MAXMSG];  // error output buffer
  
  do { // Beginning of outer ERROR handling loop.
    
    // Allocate memory:
    psc = (PointSourceCatalog*)malloc(sizeof(PointSourceCatalog));
    if (NULL==psc) {
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory available to store the PointSourceCatalog!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }
    psc->nsources = 0;
    
    psc->sources = (PointSource*)malloc(MAX_N_POINTSOURCES*sizeof(PointSource));
    if(NULL==psc->sources) {
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory available to store the PointSourceCcatalog!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    } 
    // END of memory allocation


    // At the moment the selected catalog is empty:
    PointSource ps;
    int file_counter;
    for(file_counter=0; (file_counter<psfc->nfiles)&&(EXIT_SUCCESS==*status); 
	file_counter++) {

      do {  // Beginning of inner ERROR handling loop.

	// Read-in all sources from the individual table rows, one after another:
	int source_counter; // counter to access the individual sources in the catalog
	for (source_counter=0; 
	     (source_counter < psfc->files[file_counter]->nrows) && (EXIT_SUCCESS==*status); 
	     source_counter++) {
	  // Read source data (right asension, declination, photon rate, spectrum):
	  //if (get_srctbl_row(psf->files[file_counter], source_counter, 
	  //		     psf->columns[file_counter], 
	  //		     &ra, &dec, &countrate, status))
	  if (get_PointSourceTable_Row(psfc->files[file_counter], source_counter,
				       &ps, status)) break;

	  // Get a unit vector pointing in the direction of the source:
	  Vector source_direction = unit_vector(ps.ra, ps.dec);
	  
	  // Check whether the source should be added to the preselected catalog:
	  if(fabs(scalar_product(&source_direction, &normal_vector)) < max_align) {
	    if(psc->nsources > MAX_N_POINTSOURCES) {
	      // Too many sources !
	      *status=EXIT_FAILURE;
	      sprintf(msg, "Error: too many sources (%ld)!\n", psc->nsources+1);
	      HD_ERROR_THROW(msg, *status);
	      break;
	    }

	    // Add the current source to the selected catalog:
	    psc->sources[psc->nsources].rate = ps.rate;

	    // save the source direction in the source catalog-array:
	    psc->sources[psc->nsources].ra  = ps.ra;
	    psc->sources[psc->nsources].dec = ps.dec;
	    // set lightcurve pointer to NULL
	    psc->sources[psc->nsources].lightcurve = NULL;
	    // so far there was no photon created for this source
	    psc->sources[psc->nsources].t_last_photon = -1.;
	    // source spectrum
	    psc->sources[psc->nsources].spectrum = &(spectrum_store.spectrum[0]);

	    // increase number of sources in the selected catalog
	    psc->nsources++;
	  }
	} // END of loop over all sources in the current catalog file
      } while (0); // END of inner error handling loop
    } // END of scanning the source catalogs
  } while (0); // END of outer error handling loop

  return(psc);
}



////////////////////////////////////////////////////
void free_PointSourceCatalog(PointSourceCatalog* psc)
{
  if (NULL != psc) {
    int count;
    for (count=0; count<psc->nsources; count++) {
      if (psc->sources[count].lightcurve != NULL) {
	free(psc->sources[count].lightcurve);
      }
    }

    free(psc);
  }
}



//////////////////////////////////////////////////////
int get_PointSourceTable_Row(PointSourceFile* psf, long row, PointSource* ps, int* status) 
{
  int anynul;

  // Set default values.
  ps->ra = 0.;
  ps->dec = 0.;
  ps->rate = 0.;
  ps->spectrum = NULL;
  ps->lightcurve = NULL;
  ps->t_last_photon = -1;

  // Read the data from the FITS table.

  // int fits_read_col / ffgcv
  //    (fitsfile *fptr, int datatype, int colnum, LONGLONG firstrow, LONGLONG firstelem,
  //     LONGLONG nelements, DTYPE *nulval, DTYPE *array, int *anynul, int *status) 

  // Right Ascension:
  fits_read_col(psf->fptr, TFLOAT, psf->cra, row+1, 1, 1, &ps->ra, &ps->ra, 
		&anynul, status);
  ps->ra *= M_PI/180.; // rescale from [deg]->[rad]

  // Declination:
  fits_read_col(psf->fptr, TFLOAT, psf->cdec, row+1, 1, 1, &ps->dec, &ps->dec, 
		&anynul, status);
  ps->dec *= M_PI/180.; // rescale from [deg]->[rad]

  // Photon Rate:
  fits_read_col(psf->fptr, TFLOAT, psf->crate, row+1, 1, 1, &ps->rate, &ps->rate, 
		&anynul, status);

  // TODO: Spectrum

  return(anynul);  
}



