#include "pointsources.h"


PointSourceFile* get_PointSourceFile() {
  PointSourceFile* psf = (PointSourceFile*)malloc(sizeof(PointSourceFile));

  if (NULL!=psf) {
    psf->fptr = NULL;
    psf->cra = 0;
    psf->cdec = 0;
    psf->crate = 0;
    psf->cspectrum = 0;
    psf->clightcurve = 0;
  }

  return(psf);
}



PointSourceFile* get_PointSourceFile_fromFile(char* filename, int hdu, int* status) 
{
  PointSourceFile* psf = NULL;
  char msg[MAXMSG];      // error output buffer

  do { // Beginning of ERROR handling loop.
    
    // Call the generic constructor.
    psf=get_PointSourceFile();
    if (NULL==psf) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for PointSourceFile failed!\n", *status);
      break;
    }

    // OPEN the specified FITS file and store basic information.
    headas_chat(5, "open point-source catalog in file '%s' ...\n", filename);

    // Open the source catalogue (FITS-file):
    if(fits_open_file(&psf->fptr, filename, READONLY, status)) break;
    int hdutype; // Type of the HDU
    if(fits_movabs_hdu(psf->fptr, hdu, &hdutype, status)) break;
    // Image HDU results in an error message:
    if (IMAGE_HDU==hdutype) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: FITS extension in source catalog file '%s' is "
	      "not a table but an image (HDU number: %d)!\n", filename, hdu);
      HD_ERROR_THROW(msg, *status);
      break;
    }

    // Determine the column numbers of the right ascension, declination,
    // photon rate, and spectrum columns in the FITS table
    if (fits_get_colnum(psf->fptr, CASEINSEN, "RA", &psf->cra, status)) break;
    if (fits_get_colnum(psf->fptr, CASEINSEN, "DEC", &psf->cdec, status)) break;
    if (fits_get_colnum(psf->fptr, CASEINSEN, "PPS", &psf->crate, status)) break;
    if (fits_get_colnum(psf->fptr, CASEINSEN, "SPECTRUM", &psf->cspectrum, status)) break;
    if (fits_get_colnum(psf->fptr, CASEINSEN, "LIGHTCUR", &psf->clightcurve, status)) break;

    // Determine the number of rows in the FITS table:
    if (fits_get_num_rows(psf->fptr, &psf->nrows, status)) break;	
    headas_chat(5, " contains %ld sources\n", psf->nrows);

    // Load spectra specified in the FITS header.
    *status = loadSpectra(psf->fptr, &psf->spectrumstore);
    if (EXIT_SUCCESS!=*status) break;

  } while (0); // END of ERROR handling loop.

  return(psf);
}



PointSourceFile* get_PointSourceFile_fromHDU(fitsfile* fptr, int* status) 
{
  PointSourceFile* psf = NULL;

  do { // Beginning of ERROR handling loop.
    
    // Call the generic constructor.
    psf=get_PointSourceFile();
    if (NULL==psf) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for PointSourceFile failed!\n", *status);
      break;
    }
    psf->fptr = fptr;

    // Determine the column numbers of the right ascension, declination,
    // photon rate, and spectrum columns in the FITS table
    if (fits_get_colnum(psf->fptr, CASEINSEN, "RA", &psf->cra, status)) break;
    if (fits_get_colnum(psf->fptr, CASEINSEN, "DEC", &psf->cdec, status)) break;
    if (fits_get_colnum(psf->fptr, CASEINSEN, "PPS", &psf->crate, status)) break;
    if (fits_get_colnum(psf->fptr, CASEINSEN, "SPECTRUM", &psf->cspectrum, status)) break;
    if (fits_get_colnum(psf->fptr, CASEINSEN, "LIGHTCUR", &psf->clightcurve, status)) break;

    // Determine the number of rows in the FITS table:
    if (fits_get_num_rows(psf->fptr, &psf->nrows, status)) break;	
    headas_chat(5, " contains %ld sources\n", psf->nrows);

    // Load spectra specified in the FITS header.
    *status = loadSpectra(psf->fptr, &psf->spectrumstore);
    if (EXIT_SUCCESS!=*status) break;

  } while (0); // END of ERROR handling loop.

  return(psf);
}



void free_PointSourceFile(PointSourceFile* psf) {
  if (NULL!=psf) {

    // Close the FITS file if still open.
    if (NULL!=psf->fptr) {
      int status=EXIT_SUCCESS;
      fits_close_file(psf->fptr, &status);
    }

    // Release the SpectrumStore.
    freeSpectrumStore(&psf->spectrumstore);

    free(psf);
  }
}



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



int get_PointSourceTable_Row(PointSourceFile* psf, long row, PointSource* ps, int* status) 
{
  int anynul;

  // Set default values.
  ps->ra = 0.;
  ps->dec = 0.;

  ps->rate = 0.;
  ps->lc_type = T_LC_CONSTANT;
  ps->lc = NULL;
  ps->t_last_photon = 0.;

  ps->spectrum_index = 0;
  ps->spectrum = NULL;

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

  // Spectrum type:
  fits_read_col(psf->fptr, TLONG, psf->cspectrum, row+1, 1, 1, &ps->spectrum_index, 
		&ps->spectrum_index, &anynul, status);

  // Light curve type:
  fits_read_col(psf->fptr, TLONG, psf->clightcurve, row+1, 1, 1, &ps->lc_type, 
		&ps->lc_type, &anynul, status);

  return(anynul);  
}




