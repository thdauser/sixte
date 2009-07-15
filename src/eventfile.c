#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include "eventfile.h"


int openEventlistFile(EventlistFile* ef, char* filename, int access_mode)
{
  char msg[MAXMSG];
  int status = EXIT_SUCCESS;

  headas_chat(5, "open event list file '%s' ...\n", filename);

  // Open the FITS file table for reading:
  if (fits_open_table(&ef->fptr, filename, access_mode, &status)) return(status);;

  // Get the HDU type
  int hdutype;
  if (fits_get_hdu_type(ef->fptr, &hdutype, &status)) return(status);;
  // Image HDU results in an error message.
  if (IMAGE_HDU==hdutype) {
    status=EXIT_FAILURE;
    sprintf(msg, "Error: no table extension available in FITS file '%s'!\n", filename);
    HD_ERROR_THROW(msg, status);
    return(status);
  }

  // Determine the number of rows in the event list.
  if (fits_get_num_rows(ef->fptr, &ef->nrows, &status)) return(status);
  // Set internal row counter to first row (starting at 0).
  ef->row = 0;

  return(status);
}



int closeEventlistFile(EventlistFile* ef) 
{
  int status = EXIT_SUCCESS;

  if (NULL!=ef->fptr) {
    if (fits_close_file(ef->fptr, &status)) return(status);
    ef->fptr = NULL;
  }

  return(status);
}




int EventlistFileEOF(EventlistFile* ef) {
  if (ef->row >= ef->nrows) {
    return(1);
  } else {
    return(0);
  }
}



/*
///////////////////////////////////////////////////////////////////////
// This function reads a row of data from the event list FITS table.
int get_eventlist_row(struct Eventlist_File ef, struct Event* event, int *status) 
{
  int anynul = 0;

  // time (1st column)
  event->time = 0.;
  if (0<ef.ctime)
    fits_read_col(ef.fptr, TDOUBLE, ef.ctime, ef.row+1, 1, 1, 
		  &event->time, &event->time, &anynul, status);

  // PHA channel (2nd column)
  event->pha = 0;
  if (0<ef.cpha)
    fits_read_col(ef.fptr, TLONG, ef.cpha, ef.row+1, 1, 1, 
		  &event->pha, &event->pha, &anynul, status);

  // xi (3rd column) (RAWX)
  event->xi = 0;
  if (0<ef.crawx)
    fits_read_col(ef.fptr, TINT, ef.crawx, ef.row+1, 1, 1, 
		  &event->xi, &event->xi, &anynul, status);
  event->xi -= ef.PixelOffset;

  // yi (4th column) (RAWY)
  event->yi = 0;
  if (0<ef.crawy)
    fits_read_col(ef.fptr, TINT, ef.crawy, ef.row+1, 1, 1, 
		  &event->yi, &event->yi, &anynul, status);
  event->yi -= ef.PixelOffset;

  // frame
  event->frame = 0;
  if (0<ef.cframe) 
    fits_read_col(ef.fptr, TLONG, ef.cframe, ef.row+1, 1, 1, 
		  &event->frame, &event->frame, &anynul, status);

  // patnum
  event->patnum = 0;
  if (0<ef.cpatnum)
    fits_read_col(ef.fptr, TLONG, ef.cpatnum, ef.row+1, 1, 1, 
		  &event->patnum, &event->patnum, &anynul, status);

  // patid
  event->patid = 0;
  if (0<ef.cpatid)
    fits_read_col(ef.fptr, TLONG, ef.cpatid, ef.row+1, 1, 1, 
		  &event->patid, &event->patid, &anynul, status);

  return(anynul);
}

*/
