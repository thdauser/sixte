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
    headas_chat(5, "closed event file (containing %ld rows).\n", plf->nrows);
  }

  return(status);
}


