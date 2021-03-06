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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "ladimpactfile.h"


LADImpactFile* newLADImpactFile(int* const status)
{
  LADImpactFile* file=(LADImpactFile*)malloc(sizeof(LADImpactFile));
  if (NULL==file) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for LADImpactFile failed");
    return(file);
  }

  // Initialize pointers with NULL.
  file->fptr=NULL;

  // Initialize values.
  file->nrows=0;
  file->row  =0;
  file->ctime=0;
  file->cenergy =0;
  file->cpanel  =0;
  file->cmodule =0;
  file->celement=0;
  file->cx      =0;
  file->cy      =0;
  file->cph_id  =0;
  file->csrc_id =0;

  return(file);
}


void freeLADImpactFile(LADImpactFile** const file, int* const status)
{
  if (NULL!=*file) {
    if (NULL!=(*file)->fptr) {
      fits_close_file((*file)->fptr, status);
      headas_chat(5, "closed impact list file (containing %ld rows).\n",
		  (*file)->nrows);
    }
    free(*file);
    *file=NULL;
  }
}


LADImpactFile* openLADImpactFile(const char* const filename,
				 const int mode, int* const status)
{
  LADImpactFile* file=newLADImpactFile(status);
  CHECK_STATUS_RET(*status, file);

  headas_chat(5, "open impact list file '%s' ...\n", filename);

  // Open the FITS file table for reading:
  fits_open_table(&file->fptr, filename, mode, status);
  CHECK_STATUS_RET(*status, file);

  // Get the HDU type.
  int hdutype;
  if (fits_get_hdu_type(file->fptr, &hdutype, status)) return(file);
  // Image HDU results in an error message.
  if (IMAGE_HDU==hdutype) {
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "no table extension available in FITS file '%s'",
	    filename);
    SIXT_ERROR(msg);
    return(file);
  }

  // Determine the number of rows in the impact list.
  if (fits_get_num_rows(file->fptr, &file->nrows, status)) return(file);
  // Set internal row counter to the beginning of the file (starting at 0).
  file->row = 0;

  // Determine the individual column numbers.
  fits_get_colnum(file->fptr, CASEINSEN, "TIME", &file->ctime, status);
  fits_get_colnum(file->fptr, CASEINSEN, "ENERGY", &file->cenergy, status);
  fits_get_colnum(file->fptr, CASEINSEN, "PANEL", &file->cpanel, status);
  fits_get_colnum(file->fptr, CASEINSEN, "MODULE", &file->cmodule, status);
  fits_get_colnum(file->fptr, CASEINSEN, "ELEMENT", &file->celement, status);
  fits_get_colnum(file->fptr, CASEINSEN, "X", &file->cx, status);
  fits_get_colnum(file->fptr, CASEINSEN, "Y", &file->cy, status);
  fits_get_colnum(file->fptr, CASEINSEN, "PH_ID", &file->cph_id, status);
  fits_get_colnum(file->fptr, CASEINSEN, "SRC_ID", &file->csrc_id, status);
  CHECK_STATUS_RET(*status, file);

  return(file);
}


LADImpactFile* openNewLADImpactFile(const char* const filename,
				    const char clobber,
				    int* const status)
{
  LADImpactFile* file = newLADImpactFile(status);
  CHECK_STATUS_RET(*status, file);

  // Check if the file already exists.
  int exists;
  fits_file_exists(filename, &exists, status);
  CHECK_STATUS_RET(*status, file);
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
      return(file);
    }
  }

  // Create a new event list FITS file from the template file.
  char buffer[MAXFILENAME];
  sprintf(buffer, "%s(%s%s)", filename, SIXT_DATA_PATH,
	  "/templates/ladimpactfile.tpl");
  fits_create_file(&file->fptr, buffer, status);
  CHECK_STATUS_RET(*status, file);

  // Set the time-keyword in the header.
  char datestr[MAXMSG];
  int timeref;
  fits_get_system_time(datestr, &timeref, status);
  CHECK_STATUS_RET(*status, file);
  fits_update_key(file->fptr, TSTRING, "DATE", datestr,
		  "File creation date", status);
  CHECK_STATUS_RET(*status, file);

  // Add header information about program parameters.
  // The second parameter "1" means that the headers are written
  // to the first extension.
  HDpar_stamp(file->fptr, 1, status);
  CHECK_STATUS_RET(*status, file);

  // Move to the binary table extension.
  fits_movabs_hdu(file->fptr, 2, 0, status);
  CHECK_STATUS_RET(*status, file);

  // Close the new LADImpactFile.
  freeLADImpactFile(&file, status);
  CHECK_STATUS_RET(*status, file);

  // Re-open the file.
  file=openLADImpactFile(filename, READWRITE, status);
  CHECK_STATUS_RET(*status, file);

  return(file);
}


void getNextLADImpactFromFile(LADImpactFile* const file,
			      LADImpact* const impact,
			      int* const status)
{
  // Check if the file has been opened.
  if (NULL==file) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("no impact list file opened");
    return;
  }
  if (NULL==file->fptr) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("no impact list file opened");
    return;
  }

  // Move counter to next line.
  file->row++;

  // Check if there is still a row available.
  if (file->row > file->nrows) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("impact list file contains no further entries");
    return;
  }

  // Read in the data.
  int anynul=0;
  impact->time=0.;
  fits_read_col(file->fptr, TDOUBLE, file->ctime, file->row, 1, 1,
		&impact->time, &impact->time, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->energy=0.;
  fits_read_col(file->fptr, TFLOAT, file->cenergy, file->row, 1, 1,
		&impact->energy, &impact->energy, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->panel=0;
  fits_read_col(file->fptr, TLONG, file->cpanel, file->row, 1, 1,
		&impact->panel, &impact->panel, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->module=0;
  fits_read_col(file->fptr, TLONG, file->cmodule, file->row, 1, 1,
		&impact->module, &impact->module, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->element=0;
  fits_read_col(file->fptr, TLONG, file->celement, file->row, 1, 1,
		&impact->element, &impact->element, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->position.x=0.;
  fits_read_col(file->fptr, TDOUBLE, file->cx, file->row, 1, 1,
		&impact->position.x, &impact->position.x, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->position.y=0.;
  fits_read_col(file->fptr, TDOUBLE, file->cy, file->row, 1, 1,
		&impact->position.y, &impact->position.y, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->ph_id=0;
  fits_read_col(file->fptr, TLONG, file->cph_id, file->row, 1, 1,
		&impact->ph_id, &impact->ph_id, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->src_id=0;
  fits_read_col(file->fptr, TLONG, file->csrc_id, file->row, 1, 1,
		&impact->src_id, &impact->src_id, &anynul, status);
  CHECK_STATUS_VOID(*status);

  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("reading from impact list file failed");
    return;
  }

  return;
}


void addLADImpact2File(LADImpactFile* const ilf,
		       LADImpact* const impact,
		       int* const status)
{
  ilf->row++;
  ilf->nrows++;

  fits_write_col(ilf->fptr, TDOUBLE, ilf->ctime,
		 ilf->row, 1, 1, &impact->time, status);
  fits_write_col(ilf->fptr, TFLOAT, ilf->cenergy,
		 ilf->row, 1, 1, &impact->energy, status);
  fits_write_col(ilf->fptr, TLONG, ilf->cpanel,
		 ilf->row, 1, 1, &impact->panel, status);
  fits_write_col(ilf->fptr, TLONG, ilf->cmodule,
		 ilf->row, 1, 1, &impact->module, status);
  fits_write_col(ilf->fptr, TLONG, ilf->celement,
		 ilf->row, 1, 1, &impact->element, status);
  fits_write_col(ilf->fptr, TDOUBLE, ilf->cx,
		 ilf->row, 1, 1, &(impact->position.x), status);
  fits_write_col(ilf->fptr, TDOUBLE, ilf->cy,
		 ilf->row, 1, 1, &(impact->position.y), status);
  fits_write_col(ilf->fptr, TLONG, ilf->cph_id,
		 ilf->row, 1, 1, &impact->ph_id, status);
  fits_write_col(ilf->fptr, TLONG, ilf->csrc_id,
		 ilf->row, 1, 1, &impact->src_id, status);
}
