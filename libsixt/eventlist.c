/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
*/

#include "eventlist.h"
#include "sixt.h"

int openEventList(EventList* ef, char* filename, int access_mode)
{
  char msg[MAXMSG];
  int status=EXIT_SUCCESS;

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


int closeEventList(EventList* ef) 
{
  int status=EXIT_SUCCESS;

  if (NULL!=ef) {
    if (NULL!=ef->fptr) {
      if (fits_close_file(ef->fptr, &status)) return(status);
      ef->fptr = NULL;
      headas_chat(5, "closed event file (containing %ld rows).\n", ef->nrows);
    }
  }

  return(status);
}


int EventListEOF(EventList* ef) {
  if (ef->row >= ef->nrows) {
    return(1);
  } else {
    return(0);
  }
}


int EventListRowIsValid(EventList* ef, long row) {
  // Check if the specified row is valid.
  if ((row <= 0) || (row > ef->nrows)) {
    return(0);
  } else {
    return(1);
  }
}

