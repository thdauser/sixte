#include "photonlistfile.h"


PhotonListFile* newPhotonListFile(int* const status)
{
  PhotonListFile* plf=(PhotonListFile*)malloc(sizeof(PhotonListFile));
  CHECK_NULL_RET(plf, *status, "memory allocation for PhotonListFile failed",
		 plf);

  // Initialize pointers with NULL.
  plf->fptr=NULL;

  // Initialize values.
  plf->nrows=0;
  plf->row=0;
  plf->ctime=0;
  plf->cenergy=0;
  plf->cra=0;
  plf->cdec=0;
  plf->cph_id=0;
  plf->csrc_id=0;

  return(plf);
}


void freePhotonListFile(PhotonListFile** const plf, int* const status) 
{
  if (NULL!=*plf) {
    if (NULL!=(*plf)->fptr) {
      // If the file was opened in READWRITE mode, calculate
      // the check sum an append it to the FITS header.
      int mode;
      fits_file_mode((*plf)->fptr, &mode, status);
      if (READWRITE==mode) {
	fits_write_chksum((*plf)->fptr, status);
      }
      fits_close_file((*plf)->fptr, status);
      headas_chat(5, "closed photon list file (containing %ld rows).\n", 
		  (*plf)->nrows);
    }
    free(*plf);
    *plf=NULL;
  }
}


PhotonListFile* openPhotonListFile(const char* const filename, 
				   const int access_mode,
				   int* const status)
{
  PhotonListFile* plf=newPhotonListFile(status);
  CHECK_STATUS_RET(*status, plf);
  
  headas_chat(5, "open photon list file '%s' ...\n", filename);

  // Open the FITS file table for reading:
  fits_open_table(&plf->fptr, filename, access_mode, status);
  CHECK_STATUS_RET(*status, plf);

  // Get the HDU type
  int hdutype;
  fits_get_hdu_type(plf->fptr, &hdutype, status);
  CHECK_STATUS_RET(*status, plf);
  // Image HDU results in an error message.
  if (IMAGE_HDU==hdutype) {
    char msg[MAXMSG];
    *status=EXIT_FAILURE;
    sprintf(msg, "no table extension available in file '%s'", filename);
    SIXT_ERROR(msg);
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

  // Optional columns:
  fits_write_errmark();
  int opt_status=EXIT_SUCCESS;
  fits_get_colnum(plf->fptr, CASEINSEN, "PH_ID", &plf->cph_id, &opt_status);
  opt_status=EXIT_SUCCESS;
  fits_get_colnum(plf->fptr, CASEINSEN, "SRC_ID", &plf->csrc_id, &opt_status);
  fits_clear_errmark();

  return(plf);
}


PhotonListFile* openNewPhotonListFile(const char* const filename, 
				      char* const telescop,
				      char* const instrume,
				      char* const filter,
				      const double mjdref,
				      const double timezero,
				      const double tstart,
				      const double tstop,
				      const char clobber,
				      int* const status)
{
  PhotonListFile* plf=newPhotonListFile(status);
  CHECK_STATUS_RET(*status, plf);

  // Check if the file already exists.
  int exists;
  fits_file_exists(filename, &exists, status);
  CHECK_STATUS_RET(*status, plf);
  if (0!=exists) {
    if (0!=clobber) {
      // Delete the file.
      remove(filename);
    } else {
      // Throw an error.
      char msg[MAXMSG];
      sprintf(msg, "file '%s' already exists", filename);
      SIXT_ERROR(msg);
      *status=EXIT_FAILURE;
      return(plf);
    }
  }

  // Create a new photon list FITS file from the given FITS template.
  char buffer[MAXFILENAME];
  sprintf(buffer, "%s(%s%s)", filename, SIXT_DATA_PATH, 
	  "/templates/photonlist.tpl");
  fits_create_file(&plf->fptr, buffer, status);
  CHECK_STATUS_RET(*status, plf);

  // Insert header keywords to 1st and 2nd HDU.
  sixt_add_fits_stdkeywords(plf->fptr, 1, telescop, instrume, filter,
			    mjdref, timezero, tstart, tstop, status);
  CHECK_STATUS_RET(*status, plf);
  sixt_add_fits_stdkeywords(plf->fptr, 2, telescop, instrume, filter,
			    mjdref, timezero, tstart, tstop, status);
  CHECK_STATUS_RET(*status, plf);

  // Move to the binary table extension.
  fits_movabs_hdu(plf->fptr, 2, 0, status);
  CHECK_STATUS_RET(*status, plf);

  // Close the file (it is reopened in the next step).
  freePhotonListFile(&plf, status);
  CHECK_STATUS_RET(*status, plf);

  // Open the newly created FITS file.
  plf=openPhotonListFile(filename, READWRITE, status);
  CHECK_STATUS_RET(*status, plf);

  return(plf);
}


int PhotonListFile_getNextRow(PhotonListFile* const plf, Photon* const ph)
{
  int status=EXIT_SUCCESS;

  // Move counter to next line.
  plf->row++;

  // Check if there is still a row available.
  if (plf->row > plf->nrows) {
    status=EXIT_FAILURE;
    SIXT_ERROR("photon list file contains no further entries");
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
  int anynul=0;

  // Check if there is still a row available.
  if (row > plf->nrows) {
    SIXT_ERROR("photon list file does not contain the requested line");
    return(EXIT_FAILURE);
  }

  // Read in the data.
  ph->time=0.;
  fits_read_col(plf->fptr, TDOUBLE, plf->ctime, row, 1, 1, 
		    &ph->time, &ph->time, &anynul, &status);

  ph->energy=0.;
  fits_read_col(plf->fptr, TFLOAT, plf->cenergy, row, 1, 1, 
		&ph->energy, &ph->energy, &anynul, &status);

  ph->ra=0.;
  fits_read_col(plf->fptr, TDOUBLE, plf->cra, row, 1, 1, 
		&ph->ra, &ph->ra, &anynul, &status);
  ph->ra *=M_PI/180.; // Convert from [deg] to [rad].

  ph->dec=0.;
  fits_read_col(plf->fptr, TDOUBLE, plf->cdec, row, 1, 1, 
		&ph->dec, &ph->dec, &anynul, &status);
  ph->dec*=M_PI/180.; // Convert from [deg] to [rad].

  // Optional columns: read values only if column exists.
  ph->ph_id=0;
  if (0!=plf->cph_id) {
    fits_read_col(plf->fptr, TLONG, plf->cph_id, row, 1, 1, 
		  &ph->ph_id, &ph->ph_id, &anynul, &status);
  }

  ph->src_id=0;
  if (0!=plf->csrc_id) {
    fits_read_col(plf->fptr, TLONG, plf->csrc_id, row, 1, 1, 
		  &ph->src_id, &ph->src_id, &anynul, &status);
  }

  CHECK_STATUS_RET(status, status);

  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    SIXT_ERROR("reading from event list failed");
    return(EXIT_FAILURE);
  }

  return(status);
}


int addPhoton2File(PhotonListFile* const plf, Photon* const ph)
{
  int status=EXIT_SUCCESS;

  plf->row++;
  plf->nrows++;

  // Convert from [rad] -> [deg]:
  double ra =ph->ra *180./M_PI;
  double dec=ph->dec*180./M_PI;

  // Set the unique photon identifier to the row number in the photon
  // list file. Unless the value is already set.
  if (0==ph->ph_id) {
    ph->ph_id=plf->row;
  }
  
  // Store the data in the FITS file.
  fits_write_col(plf->fptr, TDOUBLE, plf->ctime, 
		 plf->row, 1, 1, &ph->time, &status);
  fits_write_col(plf->fptr, TFLOAT, plf->cenergy, 
		 plf->row, 1, 1, &ph->energy, &status);
  fits_write_col(plf->fptr, TDOUBLE, plf->cra, 
		 plf->row, 1, 1, &ra, &status);
  fits_write_col(plf->fptr, TDOUBLE, plf->cdec, 
		 plf->row, 1, 1, &dec, &status);

  // Optional columns: write values only if column exists.
  if (0!=plf->cph_id) {
    fits_write_col(plf->fptr, TLONG, plf->cph_id, 
		   plf->row, 1, 1, &ph->ph_id, &status);
  }
  if (0!=plf->csrc_id) {
    fits_write_col(plf->fptr, TLONG, plf->csrc_id, 
		   plf->row, 1, 1, &ph->src_id, &status);
  }
  CHECK_STATUS_RET(status, status);

  return(status);
}


