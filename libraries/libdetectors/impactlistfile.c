#include "impactlistfile.h"


int openImpactListFile(ImpactListFile* ilf, char* filename, int access_mode)
{
  char msg[MAXMSG];
  int status = EXIT_SUCCESS;

  headas_chat(5, "open impact list file '%s' ...\n", filename);

  // Open the FITS file table for reading:
  if (fits_open_table(&ilf->fptr, filename, access_mode, &status)) return(status);;

  // Get the HDU type.
  int hdutype;
  if (fits_get_hdu_type(ilf->fptr, &hdutype, &status)) return(status);;
  // Image HDU results in an error message.
  if (IMAGE_HDU==hdutype) {
    status=EXIT_FAILURE;
    sprintf(msg, "Error: no table extension available in FITS file '%s'!\n", filename);
    HD_ERROR_THROW(msg, status);
    return(status);
  }

  // Determine the number of rows in the impact list.
  if (fits_get_num_rows(ilf->fptr, &ilf->nrows, &status)) return(status);
  // Set internal row counter to the beginning of the file (starting at 0).
  ilf->row = 0;


  // Determine the individual column numbers:
  // REQUIRED columns:
  if(fits_get_colnum(ilf->fptr, CASEINSEN, "TIME", &ilf->ctime, &status)) 
    return(status);
  if(fits_get_colnum(ilf->fptr, CASEINSEN, "ENERGY", &ilf->cenergy, &status)) 
    return(status);
  if(fits_get_colnum(ilf->fptr, CASEINSEN, "X", &ilf->cx, &status)) 
    return(status);
  if(fits_get_colnum(ilf->fptr, CASEINSEN, "Y", &ilf->cy, &status)) 
    return(status);

  return(status);
}



int openNewImpactListFile(ImpactListFile* ilf, char* filename, char* template)
{
  int status=EXIT_SUCCESS;

  // Set the FITS file pointer to NULL. In case that an error occurs during the file
  // generation, we want to avoid that the file pointer points somewhere.
  ilf->fptr = NULL;

  // Remove old file if it exists.
  remove(filename);

  // Create a new impact list FITS file from the given FITS template.
  fitsfile* fptr=NULL;
  char buffer[MAXMSG];
  sprintf(buffer, "%s(%s)", filename, template);
  if (fits_create_file(&fptr, buffer, &status)) return(status);

  // Close the file (it is reopened in the next step).
  if (fits_close_file(fptr, &status)) return(status);

  // Open the newly created FITS file.
  status = openImpactListFile(ilf, filename, READWRITE);

  return(status);
}



int closeImpactListFile(ImpactListFile* ilf) 
{
  int status = EXIT_SUCCESS;

  if (NULL!=ilf->fptr) {
    if (fits_close_file(ilf->fptr, &status)) return(status);
    ilf->fptr = NULL;
    headas_chat(5, "closed event file (containing %ld rows).\n", ilf->nrows);
  }

  return(status);
}



int getNextImpactListFileRow(ImpactListFile* ilf, Impact* impact) 
{
  int status = EXIT_SUCCESS;
  int anynul = 0;

  // Move counter to next line.
  ilf->row++;

  // Check if there is still a row available.
  if (ilf->row > ilf->nrows) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: impact list file contains no further entries!\n", status);
    return(status);
  }

  // Read in the data.
  impact->time = 0.;
  if (0<ilf->ctime) fits_read_col(ilf->fptr, TDOUBLE, ilf->ctime, ilf->row, 1, 1, 
				  &impact->time, &impact->time, &anynul, &status);
  impact->energy = 0.;
  if (0<ilf->cenergy) fits_read_col(ilf->fptr, TFLOAT, ilf->cenergy, ilf->row, 1, 1, 
				    &impact->energy, &impact->energy, &anynul, &status);
  impact->position.x = 0.;
  if (0<ilf->cx) fits_read_col(ilf->fptr, TDOUBLE, ilf->cx, ilf->row, 1, 1, 
			       &impact->position.x, &impact->position.x, &anynul, &status);
  impact->position.y = 0.;
  if (0<ilf->cy) fits_read_col(ilf->fptr, TDOUBLE, ilf->cy, ilf->row, 1, 1, 
			       &impact->position.y, &impact->position.y, &anynul, &status);
  
  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: reading from impact list failed!\n", status);
    return(status);
  }

  return(status);
}



int ImpactListFile_EOF(ImpactListFile* ilf) {
  if (ilf->row >= ilf->nrows) {
    return(1);
  } else {
    return(0);
  }
}


