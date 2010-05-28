#include "pointsourcefile.h"


PointSourceFile* get_PointSourceFile_fromFile(char* filename, int hdu, int* status) 
{
  fitsfile* fptr;
  char msg[MAXMSG]; // Error output buffer.
  
  // OPEN the specified FITS file and store basic information.
  headas_chat(5, "open point-source catalog in file '%s' ...\n", filename);

  // Open the source catalog (FITS-file):
  if(fits_open_file(&fptr, filename, READONLY, status)) return(NULL);;
  int hdutype; // Type of the HDU
  if(fits_movabs_hdu(fptr, hdu, &hdutype, status)) return(NULL);;
  // Image HDU results in an error message:
  if (IMAGE_HDU==hdutype) {
    *status=EXIT_FAILURE;
    sprintf(msg, "Error: FITS extension in source catalog file '%s' is "
	    "not a table but an image (HDU number: %d)!\n", filename, hdu);
    HD_ERROR_THROW(msg, *status);
    return(NULL);;
  }

  // Set up the PointSourceFile structure by using the standard routine.
  PointSourceFile* psf = get_PointSourceFile_fromHDU(fptr, status);
    
  return(psf);
}



PointSourceFile* get_PointSourceFile_fromHDU(fitsfile* fptr, int* status) 
{

  // Allocate memory for the PointSourceFile data structure.
  PointSourceFile* psf = (PointSourceFile*)malloc(sizeof(PointSourceFile));
  if (NULL==psf) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for PointSourceFile!\n",
		   *status);
    return(psf);
  }

  // Set default values.
  psf->fptr = fptr;
  psf->cra = 0;
  psf->cdec = 0;
  psf->crate = 0;
  psf->cspectrum = 0;
  psf->clightcurve = 0;
  
  // Determine the number of rows in the FITS table:
  if (fits_get_num_rows(psf->fptr, &psf->nrows, status)) return(psf);	
  headas_chat(5, " contains %ld sources\n", psf->nrows);

  // Determine the column numbers of the right ascension, declination,
  // photon rate, and spectrum columns in the FITS table
  if (fits_get_colnum(psf->fptr, CASEINSEN, "RA", &psf->cra, status)) return(psf);
  if (fits_get_colnum(psf->fptr, CASEINSEN, "DEC", &psf->cdec, status)) return(psf);
  if (fits_get_colnum(psf->fptr, CASEINSEN, "PPS", &psf->crate, status)) return(psf);
  if (fits_get_colnum(psf->fptr, CASEINSEN, "SPECTRUM", &psf->cspectrum, status)) 
    return(psf);
  if (fits_get_colnum(psf->fptr, CASEINSEN, "LIGHTCUR", &psf->clightcurve, status)) 
    return(psf);


  // Load spectra specified in the FITS header.
  *status = loadSpectra(psf->fptr, &psf->spectrumstore);
  if (EXIT_SUCCESS!=*status) return(psf);

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



int get_PointSourceTable_Row(PointSourceFile* psf, long row, 
			     PointSource* ps, int* status) 
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



SourceList* readSourceListFromPointSourceFileHDU(fitsfile* fptr, 
						 long* nelements, 
						 int* status)
{
  // Use the initialization routine to set up the PointSourceFile data structure.
  PointSourceFile* psf = get_PointSourceFile_fromHDU(fptr, status);
  if (EXIT_SUCCESS!=status) return(NULL);

  // Allocate memory for the SourceList.
  SourceList* sl = (SourceList*)malloc(psf->nrows*sizeof(struct SourceListEntry));
  if (NULL==sl) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: could not allocate memory to read PointSourceFile!\n",
		   *status);
    return(sl);
  }

  // Read all X-ray sources from the file and create a SourceList from
  // them.
  PointSource ps; // Input buffer.
  long row;
  for(row=0; row<psf->nrows; row++) {
    // Read a PointSource from the file.
    get_PointSourceTable_Row(psf, row, &ps, status);
    if (EXIT_SUCCESS!=status) return(NULL);    

    // Create an entry in the SourceList.
    sl[row].location = unit_vector(ps.ra, ps.dec);
    sl[row].line     = row;
  }

  // Return the number of elements in the SourceList to the calling function.
  *nelements = psf->nrows;

  return(sl);
}

