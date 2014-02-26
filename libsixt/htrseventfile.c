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

#include "htrseventfile.h"


int openHTRSEventFile(HTRSEventFile* hef, 
		      char* const filename,
		      const int access_mode)
{
  int status=EXIT_SUCCESS;

  // Call the corresponding routine of the underlying structure.
  status=openEventList(&hef->generic, filename, access_mode);
  CHECK_STATUS_RET(status, status);

  // Determine the HTRS-specific elements of the event list.
  // Determine the individual column numbers:
  // REQUIRED columns:
  fits_get_colnum(hef->generic.fptr, CASEINSEN, "TIME", &hef->ctime, &status);
  fits_get_colnum(hef->generic.fptr, CASEINSEN, "PHA", &hef->cpha, &status); 
  fits_get_colnum(hef->generic.fptr, CASEINSEN, "ENERGY", &hef->cenergy, &status);
  fits_get_colnum(hef->generic.fptr, CASEINSEN, "PIXEL", &hef->cpixel, &status); 
  fits_get_colnum(hef->generic.fptr, CASEINSEN, "GRADE", &hef->cgrade, &status);
  fits_get_colnum(hef->generic.fptr, CASEINSEN, "X", &hef->cx, &status);
  fits_get_colnum(hef->generic.fptr, CASEINSEN, "Y", &hef->cy, &status);
  CHECK_STATUS_RET(status, status);

  return(status);
}


int openNewHTRSEventFile(HTRSEventFile* hef, 
			 char* const filename, 
			 char* const template)
{
  int status=EXIT_SUCCESS;

  // Set the FITS file pointer to NULL. In case that an error 
  // occurs during the file generation, we want to avoid that 
  // the file pointer points somewhere.
  hef->generic.fptr = NULL;

  // Remove old file if it exists.
  remove(filename);

  // Create a new event list FITS file from a FITS template.
  fitsfile* fptr=NULL;
  char buffer[MAXMSG];
  sprintf(buffer, "%s(%s)", filename, template);
  fits_create_file(&fptr, buffer, &status);
  CHECK_STATUS_RET(status, status);

  // Set the time-keyword in the event list header.
  char datestr[MAXMSG];
  int timeref;
  fits_get_system_time(datestr, &timeref, &status);
  CHECK_STATUS_RET(status, status);
  fits_update_key(fptr, TSTRING, "DATE", datestr, 
		  "File creation date", &status);
  CHECK_STATUS_RET(status, status);

  // Add header information about program parameters.
  // The second parameter "1" means that the headers are writte
  // to the first extension.
  HDpar_stamp(fptr, 1, &status);
  CHECK_STATUS_RET(status, status);

  // Close the file. It will be re-opened immediately with the
  // standard opening routine.
  fits_close_file(fptr, &status);
  CHECK_STATUS_RET(status, status);

  // Open the newly created FITS file.
  status=openHTRSEventFile(hef, filename, READWRITE);

  return(status);
}


int closeHTRSEventFile(HTRSEventFile* hef)
{
  // Call the corresponding routine of the underlying structure.
  return(closeEventList(&hef->generic));
}


int addHTRSEvent2File(HTRSEventFile* hef, HTRSEvent* event)
{
  int status=EXIT_SUCCESS;

  // Insert a new, empty row to the table:
  fits_insert_rows(hef->generic.fptr, hef->generic.row, 1, &status);
  CHECK_STATUS_RET(status, status);
  hef->generic.row++;
  hef->generic.nrows++;

  // Write the event data to the newly created row.
  status=HTRSEventFile_writeRow(hef, event, hef->generic.row);

  return(status);
}


int HTRSEventFile_getNextRow(HTRSEventFile* hef, HTRSEvent* event)
{
  int status=EXIT_SUCCESS;

  // Move counter to next line.
  hef->generic.row++;

  // Check if there is still a row available.
  if (hef->generic.row > hef->generic.nrows) {
    status=EXIT_FAILURE;
    SIXT_ERROR("event list contains no further entries");
    return(status);
  }

  // Read the new HTRSEvent from the file.
  status=HTRSEventFile_getRow(hef, event, hef->generic.row);

  return(status);
}


int HTRSEventFile_getRow(HTRSEventFile* hef, HTRSEvent* event, const long row)
{
  int status=EXIT_SUCCESS;
  int anynul=0;

  // Check if there is still a row available.
  if (row > hef->generic.nrows) {
    status=EXIT_FAILURE;
    SIXT_ERROR("event list does not contain the requested line");
    return(status);
  }

  // Read in the data.
  event->time=0.;
  fits_read_col(hef->generic.fptr, TDOUBLE, hef->ctime, row, 1, 1, 
		&event->time, &event->time, &anynul, &status);
  event->pha=0;
  fits_read_col(hef->generic.fptr, TLONG, hef->cpha, row, 1, 1, 
		&event->pha, &event->pha, &anynul, &status);
  event->energy=0;
  fits_read_col(hef->generic.fptr, TFLOAT, hef->cenergy, row, 1, 1, 
		&event->energy, &event->energy, &anynul, &status);
  event->pixel=0;
  fits_read_col(hef->generic.fptr, TINT, hef->cpixel, row, 1, 1, 
		&event->pixel, &event->pixel, &anynul, &status);
  event->pixel--;
  
  event->grade=0;
  fits_read_col(hef->generic.fptr, TINT, hef->cgrade, row, 1, 1, 
		&event->grade, &event->grade, &anynul, &status);

  event->x=0.;
  fits_read_col(hef->generic.fptr, TDOUBLE, hef->cx, row, 1, 1, 
		&event->x, &event->x, &anynul, &status);
  event->y=0.;
  fits_read_col(hef->generic.fptr, TDOUBLE, hef->cy, row, 1, 1, 
		&event->y, &event->y, &anynul, &status);
  CHECK_STATUS_RET(status, status);
  
  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    status=EXIT_FAILURE;
    SIXT_ERROR("failed reading from event list");
    return(status);
  }

  return(status);
}


int HTRSEventFile_writeRow(HTRSEventFile* hef, HTRSEvent* event, const long row) {
  int status=EXIT_SUCCESS;

  fits_write_col(hef->generic.fptr, TDOUBLE, hef->ctime, row, 
		 1, 1, &event->time, &status);
  fits_write_col(hef->generic.fptr, TLONG, hef->cpha, row, 
		 1, 1, &event->pha, &status);
  fits_write_col(hef->generic.fptr, TFLOAT, hef->cenergy, row, 
		 1, 1, &event->energy, &status);

  int pixel=event->pixel+1;
  fits_write_col(hef->generic.fptr, TINT, hef->cpixel, row, 
		 1, 1, &pixel, &status);

  fits_write_col(hef->generic.fptr, TINT, hef->cgrade, row, 
		 1, 1, &event->grade, &status);
  fits_write_col(hef->generic.fptr, TDOUBLE, hef->cx, row, 
		 1, 1, &event->x, &status);
  fits_write_col(hef->generic.fptr, TDOUBLE, hef->cy, row, 
		 1, 1, &event->y, &status);
  CHECK_STATUS_RET(status, status);

  return(status);
}

