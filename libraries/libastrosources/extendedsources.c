#include "extendedsources.h"


ExtendedSourceFile* get_ExtendedSourceFile() {
  ExtendedSourceFile* esf = 
    (ExtendedSourceFile*)malloc(sizeof(ExtendedSourceFile));

  if (NULL!=esf) {
    esf->fptr = NULL;
    esf->cra = 0;
    esf->cdec = 0;
    esf->cradius = 0;
    esf->crate = 0;
    esf->cspectrum = 0;
    esf->clightcurve = 0;
  }

  return(esf);
}



ExtendedSourceFile* get_ExtendedSourceFile_fromFile(char* filename, 
						    int hdu, 
						    int* status) 
{
  ExtendedSourceFile* esf = NULL;
  char msg[MAXMSG]; // Error output buffer.

  do { // Beginning of ERROR handling loop.
    
    // Call the generic constructor.
    esf=get_ExtendedSourceFile();
    if (NULL==esf) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for ExtendedSourceFile failed!\n", 
		     *status);
      break;
    }

    // OPEN the specified FITS file and store basic information.
    headas_chat(5, "open extended-source catalog in file '%s' ...\n", filename);

    // Open the source catalogue (FITS-file):
    if(fits_open_file(&esf->fptr, filename, READONLY, status)) break;
    int hdutype; // Type of the HDU
    if(fits_movabs_hdu(esf->fptr, hdu, &hdutype, status)) break;
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
    if (fits_get_colnum(esf->fptr, CASEINSEN, "RA", &esf->cra, status)) break;
    if (fits_get_colnum(esf->fptr, CASEINSEN, "DEC", &esf->cdec, status)) break;
    if (fits_get_colnum(esf->fptr, CASEINSEN, "RADIUS", &esf->cradius, status)) break;
    if (fits_get_colnum(esf->fptr, CASEINSEN, "PPS", &esf->crate, status)) break;
    if (fits_get_colnum(esf->fptr, CASEINSEN, "SPECTRUM", &esf->cspectrum, status)) break;
    if (fits_get_colnum(esf->fptr, CASEINSEN, "LIGHTCUR", &esf->clightcurve, status)) break;

    // Determine the number of rows in the FITS table:
    if (fits_get_num_rows(esf->fptr, &esf->nrows, status)) break;	
    headas_chat(5, " contains %ld sources\n", esf->nrows);

    // Load spectra specified in the FITS header.
    *status = loadSpectra(esf->fptr, &esf->spectrumstore);
    if (EXIT_SUCCESS!=*status) break;

  } while (0); // END of ERROR handling loop.

  return(esf);
}



ExtendedSourceFile* get_ExtendedSourceFile_fromHDU(fitsfile* fptr, int* status) 
{
  ExtendedSourceFile* esf = NULL;

  do { // Beginning of ERROR handling loop.
    
    // Call the generic constructor.
    esf=get_ExtendedSourceFile();
    if (NULL==esf) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for ExtendedSourceFile failed!\n", *status);
      break;
    }
    esf->fptr = fptr;

    // Determine the column numbers of the right ascension, declination,
    // photon rate, and spectrum columns in the FITS table
    if (fits_get_colnum(esf->fptr, CASEINSEN, "RA", &esf->cra, status)) break;
    if (fits_get_colnum(esf->fptr, CASEINSEN, "DEC", &esf->cdec, status)) break;
    if (fits_get_colnum(esf->fptr, CASEINSEN, "RADIUS", &esf->cradius, status)) break;
    if (fits_get_colnum(esf->fptr, CASEINSEN, "PPS", &esf->crate, status)) break;
    if (fits_get_colnum(esf->fptr, CASEINSEN, "SPECTRUM", &esf->cspectrum, status)) break;
    if (fits_get_colnum(esf->fptr, CASEINSEN, "LIGHTCUR", &esf->clightcurve, status)) break;

    // Determine the number of rows in the FITS table:
    if (fits_get_num_rows(esf->fptr, &esf->nrows, status)) break;	
    headas_chat(5, " contains %ld sources\n", esf->nrows);

    // Load spectra specified in the FITS header.
    *status = loadSpectra(esf->fptr, &esf->spectrumstore);
    if (EXIT_SUCCESS!=*status) break;

  } while (0); // END of ERROR handling loop.

  return(esf);
}



void free_ExtendedSourceFile(ExtendedSourceFile* esf) {
  if (NULL!=esf) {

    // Close the FITS file if still open.
    if (NULL!=esf->fptr) {
      int status=EXIT_SUCCESS;
      fits_close_file(esf->fptr, &status);
    }

    // Release the SpectrumStore.
    freeSpectrumStore(&esf->spectrumstore);

    free(esf);
  }
}



ExtendedSourceCatalog* getExtendedSourceCatalog(ExtendedSourceFile* esf, 
						Vector normal_vector, 
						const double max_align,
						int* status)
{
  ExtendedSourceCatalog* esc = NULL;
  char msg[MAXMSG];  // Error output buffer.
  
  do { // Beginning of outer ERROR handling loop.
    
    // Allocate memory:
    esc = (ExtendedSourceCatalog*)malloc(sizeof(ExtendedSourceCatalog));
    if (NULL==esc) {
      *status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: not enough memory available to store the ExtendedSourceCatalog!\n", 
		     *status);
      break;
    }
    esc->nsources = 0;
    
    esc->sources = (ExtendedSource*)malloc(MAX_N_EXTENDEDSOURCES*sizeof(ExtendedSource));
    if(NULL==esc->sources) {
      *status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: not enough memory available to store the ExtendedSourceCcatalog!\n", 
		     *status);
      break;
    } 
    // END of memory allocation


    // At the moment the selected catalog is empty.
    ExtendedSource es;
    // Counter to access the individual sources in the catalog.
    int source_counter; 
    // Unit vector pointing in the direction of the source.
    Vector source_direction;

    // Read-in all sources from the individual table rows, one after another:
    for (source_counter=0; source_counter < esf->nrows; source_counter++) {
      
      // Read source data (right asension, declination, photon rate, spectrum):
      if (get_ExtendedSourceTable_Row(esf, source_counter, &es, status)) break;

      // Get a unit vector pointing in the direction of the source:
      source_direction = unit_vector(es.ra, es.dec);
	  
      // Check whether the source should be added to the preselected catalog:
      if(fabs(scalar_product(&source_direction, &normal_vector)) < max_align) {
	if(esc->nsources > MAX_N_EXTENDEDSOURCES) {
	  // Too many sources in the ExtendedSourceCatalog !
	  *status=EXIT_FAILURE;
	  sprintf(msg, "Error: too many sources (%ld) in the ExtendedSourceCatalog!\n", 
		  esc->nsources+1);
	  HD_ERROR_THROW(msg, *status);
	  break;
	}

	// Add the current source to the selected catalog:
	esc->sources[esc->nsources].ra  = es.ra;
	esc->sources[esc->nsources].dec = es.dec;
	esc->sources[esc->nsources].radius = es.radius;
	esc->sources[esc->nsources].rate = es.rate;

	// Set the light curve type to the value specified in the source catalog.
	esc->sources[esc->nsources].lc_type = es.lc_type;
	// If the light curve type is a value > 0 meaning that the particular
	// light curve for this source is given in a FITS file, load that file
	// and assign the pointer to the light curve.
	if (0<es.lc_type) {
	  // Try and find out the name of the FITS file containing the light curve.
	  char lc_filename[MAXMSG], key[MAXMSG], comment[MAXMSG]; // Buffers
	  sprintf(key, "LC%06ld", es.lc_type);
	  if (fits_read_key(esf->fptr, TSTRING, key, 
			    lc_filename, comment, status)) break;
	  // Then load the file data and store the light curve.
	  esc->sources[esc->nsources].lc = 
	    loadLinLightCurveFromFile(lc_filename, es.rate, status);
	  if (EXIT_SUCCESS!=*status) break;
	  // So far there was no photon created for this source.
	  esc->sources[esc->nsources].t_last_photon = 0.;
	  
	} else {
	  // No particular light curve has been assigned to this source.
	  // Set light curve pointer to NULL.
	  esc->sources[esc->nsources].lc = NULL;
	  // So far there was no photon created for this source.
	  esc->sources[esc->nsources].t_last_photon = 0.;
	  
	} // END of assigning a light curve to the source.
	  
	// Source spectrum.
	if ((es.spectrum_index<1) || 
	    (es.spectrum_index>esf->spectrumstore.nspectra)) {
	  headas_chat(0, "\n### Warning: no source spectrum specified for extended source!\n"
		      "     Using default spectrum instead.\n");
	  esc->sources[esc->nsources].spectrum = &esf->spectrumstore.spectrum[0];
	} else {
	  esc->sources[esc->nsources].spectrum = 
	    &esf->spectrumstore.spectrum[es.spectrum_index-1];
	} // END of assigning a spectrum to the source.
	
	// Increase number of sources in the selected catalog
	esc->nsources++;
      } // END of check whether this sources lies within the preselection band.
    } // END of loop over all sources in the current catalog file
    if (EXIT_SUCCESS!=*status) break;
  } while (0); // END of outer error handling loop

  return(esc);
}



void free_ExtendedSourceCatalog(ExtendedSourceCatalog* esc)
{
  if (NULL!=esc) {
    if (NULL!=esc->sources) {
      int count;
      for (count=0; count<esc->nsources; count++) {
	if (esc->sources[count].lc != NULL) {
	  freeLinLightCurve(esc->sources[count].lc);
	  esc->sources[count].lc=NULL;
	}
      }
      free(esc->sources);
      esc->sources=NULL;
    }
    free(esc);
  }
}



int get_ExtendedSourceTable_Row(ExtendedSourceFile* esf, long row, ExtendedSource* es, int* status) 
{
  int anynul;

  // Set default values.
  es->ra = 0.;
  es->dec = 0.;
  es->radius = 0.;
  es->rate = 0.;
  es->lc_type = T_LC_CONSTANT;
  es->lc = NULL;
  es->t_last_photon = 0.;

  es->spectrum_index = 0;
  es->spectrum = NULL;

  // Read the data from the FITS table.

  // int fits_read_col / ffgcv
  //    (fitsfile *fptr, int datatype, int colnum, LONGLONG firstrow, LONGLONG firstelem,
  //     LONGLONG nelements, DTYPE *nulval, DTYPE *array, int *anynul, int *status) 

  // Right Ascension:
  fits_read_col(esf->fptr, TFLOAT, esf->cra, row+1, 1, 1, &es->ra, &es->ra, 
		&anynul, status);
  es->ra *= M_PI/180.; // rescale from [deg]->[rad]

  // Declination:
  fits_read_col(esf->fptr, TFLOAT, esf->cdec, row+1, 1, 1, &es->dec, &es->dec, 
		&anynul, status);
  es->dec *= M_PI/180.; // rescale from [deg]->[rad]

  // Source extension (radius):
  fits_read_col(esf->fptr, TFLOAT, esf->cradius, row+1, 1, 1, &es->radius, &es->radius, 
		&anynul, status);
  es->radius *= M_PI/180.; // rescale from [deg]->[rad]

  // Photon Rate:
  fits_read_col(esf->fptr, TFLOAT, esf->crate, row+1, 1, 1, &es->rate, &es->rate, 
		&anynul, status);

  // Spectrum type:
  fits_read_col(esf->fptr, TLONG, esf->cspectrum, row+1, 1, 1, &es->spectrum_index, 
		&es->spectrum_index, &anynul, status);

  // Light curve type:
  fits_read_col(esf->fptr, TLONG, esf->clightcurve, row+1, 1, 1, &es->lc_type, 
		&es->lc_type, &anynul, status);

  return(anynul);  
}




