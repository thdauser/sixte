#include "photonlistfile.h"


int openPhotonListFile(PhotonListFile* plf, char* filename, int access_mode)
{
  char msg[MAXMSG];
  int status = EXIT_SUCCESS;

  headas_chat(5, "open photon list file '%s' ...\n", filename);

  // Open the FITS file table for reading:
  if (fits_open_table(&plf->fptr, filename, access_mode, &status)) return(status);;

  // Get the HDU type
  int hdutype;
  if (fits_get_hdu_type(plf->fptr, &hdutype, &status)) return(status);;
  // Image HDU results in an error message.
  if (IMAGE_HDU==hdutype) {
    status=EXIT_FAILURE;
    sprintf(msg, "Error: no table extension available in FITS file '%s'!\n", filename);
    HD_ERROR_THROW(msg, status);
    return(status);
  }

  // Determine the number of rows in the photon list.
  if (fits_get_num_rows(plf->fptr, &plf->nrows, &status)) return(status);
  // Set internal row counter to the beginning of the file (starting at 0).
  plf->row = 0;


  // Determine the individual column numbers:
  // REQUIRED columns:
  if(fits_get_colnum(plf->fptr, CASEINSEN, "TIME", &plf->ctime, &status)) 
    return(status);
  if(fits_get_colnum(plf->fptr, CASEINSEN, "ENERGY", &plf->cenergy, &status)) 
    return(status);
  if(fits_get_colnum(plf->fptr, CASEINSEN, "RA", &plf->cra, &status)) 
    return(status);
  if(fits_get_colnum(plf->fptr, CASEINSEN, "DEC", &plf->cdec, &status)) 
    return(status);

  return(status);
}



int openNewPhotonListFile(PhotonListFile* plf, char* filename, char* template)
{
  int status=EXIT_SUCCESS;

  // Set the FITS file pointer to NULL. In case that an error occurs during the file
  // generation, we want to avoid that the file pointer points somewhere.
  plf->fptr = NULL;

  // Remove old file if it exists.
  remove(filename);

  // Create a new photon list FITS file from the given FITS template.
  fitsfile* fptr=NULL;
  char buffer[MAXMSG];
  sprintf(buffer, "%s(%s)", filename, template);
  if (fits_create_file(&fptr, buffer, &status)) return(status);

  // Close the file (it is reopened in the next step).
  if (fits_close_file(fptr, &status)) return(status);

  // Open the newly created FITS file.
  status = openPhotonListFile(plf, filename, READWRITE);

  return(status);
}



int closePhotonListFile(PhotonListFile* plf) 
{
  int status = EXIT_SUCCESS;

  if (NULL!=plf->fptr) {
    if (fits_close_file(plf->fptr, &status)) return(status);
    plf->fptr = NULL;
    headas_chat(5, "closed photon list file (containing %ld rows).\n", plf->nrows);
  }

  return(status);
}



int PhotonListFile_getNextRow(PhotonListFile* plf, Photon* ph)
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



int PhotonListFile_getRow(PhotonListFile* plf, Photon* ph, long row)
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

  
  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: reading from event list failed!\n", status);
    return(status);
  }

  return(status);
}



int addPhoton2File(PhotonListFile* plf, Photon* ph)
{
  int status=EXIT_SUCCESS;

  // Insert a new, empty row to the table:
  if (fits_insert_rows(plf->fptr, plf->row, 1, &status)) return(status);
  plf->row++;
  plf->nrows++;

  if (fits_write_col(plf->fptr, TDOUBLE, plf->ctime, 
		     plf->row, 1, 1, &ph->time, &status)) return(status);
  if (fits_write_col(plf->fptr, TFLOAT, plf->cenergy, 
		     plf->row, 1, 1, &ph->energy, &status)) return(status);
  if (fits_write_col(plf->fptr, TDOUBLE, plf->cra, 
		     plf->row, 1, 1, &ph->ra, &status)) return(status);
  if (fits_write_col(plf->fptr, TDOUBLE, plf->cdec, 
		     plf->row, 1, 1, &ph->dec, &status)) return(status);

  return(status);
}


