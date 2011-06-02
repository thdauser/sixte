#include "photonlistfile.h"


PhotonListFile* openPhotonListFile(const char* const filename, 
				   const int access_mode,
				   int* const status)
{
  PhotonListFile* plf=NULL;
  char msg[MAXMSG];

  headas_chat(5, "open photon list file '%s' ...\n", filename);

  // Get memory for the new PhotonListFile object.
  plf = (PhotonListFile*)malloc(sizeof(PhotonListFile));
  if (NULL==plf) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("memory allocation for PhotonListFile failed!\n", *status);
    return(NULL);
  }
  
  // Initialize with default values.
  plf->fptr=NULL;
  plf->nrows=0;
  plf->row=0;
  plf->ctime=0;
  plf->cenergy=0;
  plf->cra=0;
  plf->cdec=0;
  plf->cph_id=0;
  plf->csrc_id=0;

  // Open the FITS file table for reading:
  fits_open_table(&plf->fptr, filename, access_mode, status);
  CHECK_STATUS_RET(*status, plf);

  // Get the HDU type
  int hdutype;
  fits_get_hdu_type(plf->fptr, &hdutype, status);
  CHECK_STATUS_RET(*status, plf);
  // Image HDU results in an error message.
  if (IMAGE_HDU==hdutype) {
    *status=EXIT_FAILURE;
    sprintf(msg, "Error: no table extension available in FITS file '%s'!\n", filename);
    HD_ERROR_THROW(msg, *status);
    return(plf);
  }

  // Determine the number of rows in the photon list.
  fits_get_num_rows(plf->fptr, &plf->nrows, status);
  CHECK_STATUS_RET(*status, plf);

  // Determine the individual column numbers.
  fits_get_colnum(plf->fptr, CASEINSEN, "TIME", &plf->ctime, status); 
  CHECK_STATUS_RET(*status, plf);
  fits_get_colnum(plf->fptr, CASEINSEN, "ENERGY", &plf->cenergy, status); 
  CHECK_STATUS_RET(*status, plf);
  fits_get_colnum(plf->fptr, CASEINSEN, "RA", &plf->cra, status);
  CHECK_STATUS_RET(*status, plf);
  fits_get_colnum(plf->fptr, CASEINSEN, "DEC", &plf->cdec, status);
  CHECK_STATUS_RET(*status, plf);
  fits_get_colnum(plf->fptr, CASEINSEN, "PH_ID", &plf->cph_id, status);
  CHECK_STATUS_RET(*status, plf);
  fits_get_colnum(plf->fptr, CASEINSEN, "SRC_ID", &plf->csrc_id, status);
  CHECK_STATUS_RET(*status, plf);

  return(plf);
}


PhotonListFile* openNewPhotonListFile(const char* const filename, 
				      const char* const template,
				      int* const status)
{
  PhotonListFile* plf=NULL;

  // Remove old file, if it exists.
  remove(filename);

  // Create a new photon list FITS file from the given FITS template.
  fitsfile* fptr=NULL;
  char buffer[MAXMSG];
  sprintf(buffer, "%s(%s)", filename, template);
  fits_create_file(&fptr, buffer, status);
  CHECK_STATUS_RET(*status, plf);

  // Add header information about program parameters.
  // The second parameter "1" means that the headers are written
  // to the first extension.
  HDpar_stamp(fptr, 1, status);
  if (EXIT_SUCCESS!=*status) return(plf);

  // Move to the binary table extension.
  fits_movabs_hdu(fptr, 2, 0, status);
  CHECK_STATUS_RET(*status, plf);

  // Close the file (it is reopened in the next step).
  fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, plf);

  // Open the newly created FITS file.
  plf=openPhotonListFile(filename, READWRITE, status);
  CHECK_STATUS_RET(*status, plf);

  return(plf);
}


void freePhotonListFile(PhotonListFile** const plf, int* const status) 
{
  if (NULL!=*plf) {
    if (NULL!=(*plf)->fptr) {
      fits_close_file((*plf)->fptr, status);
      headas_chat(5, "closed photon list file (containing %ld rows).\n", 
		  (*plf)->nrows);
    }
    free(*plf);
    *plf=NULL;
  }
}


int PhotonListFile_getNextRow(PhotonListFile* const plf, Photon* const ph)
{
  int status=EXIT_SUCCESS;

  // Move counter to next line.
  plf->row++;

  // Check if there is still a row available.
  if (plf->row > plf->nrows) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: photon list file contains no further entries!\n", status);
    return(status);
  }

  // Read the new Photon from the file.
  status=PhotonListFile_getRow(plf, ph, plf->row);

  return(status);
}


int PhotonListFile_getRow(PhotonListFile* const plf, 
			  Photon* const ph, const long row)
{
  int status=EXIT_SUCCESS;
  int anynul = 0;

  // Check if there is still a row available.
  if (row > plf->nrows) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: photon list file does not contain the requested line!\n", 
		   status);
    return(status);
  }

  // Read in the data.
  ph->time = 0.;
  if (fits_read_col(plf->fptr, TDOUBLE, plf->ctime, plf->row, 1, 1, 
		    &ph->time, &ph->time, &anynul, &status)) return(status);

  ph->energy = 0.;
  if (fits_read_col(plf->fptr, TFLOAT, plf->cenergy, plf->row, 1, 1, 
		    &ph->energy, &ph->energy, &anynul, &status)) return(status);
  ph->ra = 0.;
  if (fits_read_col(plf->fptr, TDOUBLE, plf->cra, plf->row, 1, 1, 
		    &ph->ra, &ph->ra, &anynul, &status)) return(status);
  ph->dec = 0.;
  if (fits_read_col(plf->fptr, TDOUBLE, plf->cdec, plf->row, 1, 1, 
		    &ph->dec, &ph->dec, &anynul, &status)) return(status);

  ph->ph_id = 0;
  if (fits_read_col(plf->fptr, TLONG, plf->cph_id, plf->row, 1, 1, 
		    &ph->ph_id, &ph->ph_id, &anynul, &status)) return(status);
  ph->src_id = 0;
  if (fits_read_col(plf->fptr, TLONG, plf->csrc_id, plf->row, 1, 1, 
		    &ph->src_id, &ph->src_id, &anynul, &status)) return(status);

  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: reading from event list failed!\n", status);
    return(status);
  }

  // Convert from [deg] to [rad].
  ph->ra  *= M_PI/180.;
  ph->dec *= M_PI/180.;
  
  return(status);
}



int addPhoton2File(PhotonListFile* const plf, Photon* const ph)
{
  int status=EXIT_SUCCESS;

  // Convert from [rad] -> [deg]:
  double ra  = ph->ra  * 180./M_PI;
  double dec = ph->dec * 180./M_PI;

  // Insert a new, empty row to the table:
  if (fits_insert_rows(plf->fptr, plf->row, 1, &status)) return(status);
  plf->row++;
  plf->nrows++;

  // Set the unique photon identifier to the row number in the photon
  // list file. Up to now it should contain the default value 0.
  assert(0==ph->ph_id);
  ph->ph_id = plf->row;
  
  // Store the data in the FITS file.
  if (fits_write_col(plf->fptr, TDOUBLE, plf->ctime, 
		     plf->row, 1, 1, &ph->time, &status)) return(status);
  if (fits_write_col(plf->fptr, TFLOAT, plf->cenergy, 
		     plf->row, 1, 1, &ph->energy, &status)) return(status);
  if (fits_write_col(plf->fptr, TDOUBLE, plf->cra, 
		     plf->row, 1, 1, &ra, &status)) return(status);
  if (fits_write_col(plf->fptr, TDOUBLE, plf->cdec, 
		     plf->row, 1, 1, &dec, &status)) return(status);
  if (fits_write_col(plf->fptr, TLONG, plf->cph_id, 
		     plf->row, 1, 1, &ph->ph_id, &status)) return(status);
  if (fits_write_col(plf->fptr, TLONG, plf->csrc_id, 
		     plf->row, 1, 1, &ph->src_id, &status)) return(status);

  return(status);
}


