#include "impactlist.h"


int impactlist_openFile(struct ImpactlistFile* imf, char* filename, int access_mode) 
{
  char msg[MAXMSG]; // buffer for error messages
  int status = EXIT_SUCCESS;

  do {  // ERROR handling loop

    // Open the FITS file table for reading:
    if (fits_open_table(&imf->fptr, filename, access_mode, &status)) break;

    // Get the HDU type
    int hdutype;
    if (fits_get_hdu_type(imf->fptr, &hdutype, &status)) break;

    // Image HDU results in an error message.
    if (IMAGE_HDU==hdutype) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: no table extension available in impact list "
	      "FITS file '%s'!\n", filename);
      HD_ERROR_THROW(msg, status);
      break;
    }

    // Determine the number of rows in the event list.
    if (fits_get_num_rows(imf->fptr, &imf->nrows, &status)) break;

    // Set internal row counter to first row (starting at 0).
    imf->row = 0;

    // Determine the individual column numbers:
    if(fits_get_colnum(imf->fptr, CASEINSEN, "TIME", &imf->ctime, &status)) break;
    if(fits_get_colnum(imf->fptr, CASEINSEN, "ENERGY", &imf->cenergy, &status)) break;
    if(fits_get_colnum(imf->fptr, CASEINSEN, "X", &imf->cx, &status)) break;
    if(fits_get_colnum(imf->fptr, CASEINSEN, "Y", &imf->cy, &status)) break;

  } while(0);  // END of error handling loop

  return(status);
}



int impactlist_getNextRow(struct ImpactlistFile* imf, struct Impact* impact) {
  int status = EXIT_SUCCESS;
  int anynul = 0;

  // Move counter to next line.
  imf->row++;

  // Check if there is still a row available.
  if (imf->row > imf->nrows) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: impact list file contains no further entries!\n", status);
    return(status);
  }

  // Read in the data.
  // time (1st column)
  impact->time = 0.;
  if (0<imf->ctime) fits_read_col(imf->fptr, TDOUBLE, imf->ctime, imf->row, 1, 1, 
				  &impact->time, &impact->time, &anynul, &status);
  impact->energy = 0.;
  if (0<imf->cenergy) fits_read_col(imf->fptr, TFLOAT, imf->cenergy, imf->row, 1, 1, 
				    &impact->energy, &impact->energy, &anynul, &status);
  impact->position.x = 0.;
  if (0<imf->cx) fits_read_col(imf->fptr, TDOUBLE, imf->cx, imf->row, 1, 1, 
			       &impact->position.x, &impact->position.x, &anynul, &status);
  impact->position.y = 0.;
  if (0<imf->cy) fits_read_col(imf->fptr, TDOUBLE, imf->cy, imf->row, 1, 1, 
			       &impact->position.y, &impact->position.y, &anynul, &status);
  
  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: reading from impact list failed!\n", status);
    return(status);
  }

  // Read header keywords from the FITS file.
  char comment[MAXMSG];
  if ((fits_read_key(imf->fptr, TSTRING, "ATTITUDE", &imf->attitude_filename, comment, &status)))
    return(status);

  return(status);
}

