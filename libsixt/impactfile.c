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

#include "impactfile.h"


ImpactFile* newImpactFile(int* const status)
{
  ImpactFile* file=(ImpactFile*)malloc(sizeof(ImpactFile));
  if (NULL==file) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for ImpactFile failed");
    return(file);
  }

  // Initialize pointers with NULL.
  file->fptr=NULL;

  // Initialize values.
  file->nrows=0;
  file->row  =0;
  file->ctime=0;
  file->cenergy=0;
  file->cx   =0;
  file->cy   =0;
  file->cph_id =0;
  file->csrc_id=0;

  return(file);
}


void freeImpactFile(ImpactFile** const file, int* const status)
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


ImpactFile* openImpactFile(const char* const filename,
			   const int mode, int* const status)
{
  ImpactFile* file=newImpactFile(status);
  CHECK_STATUS_RET(*status, file);

  headas_chat(5, "open impact list file '%s' ...\n", filename);

  // Open the FITS file table for reading:
  if (fits_open_table(&file->fptr, filename, mode, status)) return(file);;

  // Get the HDU type.
  int hdutype;
  if (fits_get_hdu_type(file->fptr, &hdutype, status)) return(file);;
  // Image HDU results in an error message.
  if (IMAGE_HDU==hdutype) {
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "no table extension available in file '%s'", 
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
  fits_get_colnum(file->fptr, CASEINSEN, "X", &file->cx, status);
  fits_get_colnum(file->fptr, CASEINSEN, "Y", &file->cy, status);
  fits_get_colnum(file->fptr, CASEINSEN, "PH_ID", &file->cph_id, status);
  fits_get_colnum(file->fptr, CASEINSEN, "SRC_ID", &file->csrc_id, status);
  CHECK_STATUS_RET(*status, file);

  return(file);
}


ImpactFile* openNewImpactFile(const char* const filename,
			      char* const telescop,
			      char* const instrume,
			      char* const filter,
			      char* const ancrfile,
			      char* const respfile,
			      const double mjdref,
			      const double timezero,
			      const double tstart,
			      const double tstop,
			      const char clobber,
			      int* const status)
{
  ImpactFile* file=newImpactFile(status);
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

  // Create a new impact list FITS file from the template.
  char buffer[MAXFILENAME];
  sprintf(buffer, "%s(%s%s)", filename, SIXT_DATA_PATH, 
	  "/templates/impactfile.tpl");
  fits_create_file(&file->fptr, buffer, status);
  CHECK_STATUS_RET(*status, file);

  // Insert header keywords to 1st and 2nd HDU.
  sixt_add_fits_stdkeywords(file->fptr, 1, telescop, instrume, filter,
			    ancrfile, respfile,
			    mjdref, timezero, tstart, tstop, status);
  CHECK_STATUS_RET(*status, NULL);
  sixt_add_fits_stdkeywords(file->fptr, 2, telescop, instrume, filter,
			    ancrfile, respfile,
			    mjdref, timezero, tstart, tstop, status);
  CHECK_STATUS_RET(*status, NULL);

  // Move to the binary table extension.
  fits_movabs_hdu(file->fptr, 2, 0, status);
  CHECK_STATUS_RET(*status, file);

  // Close the new ImpactFile.
  freeImpactFile(&file, status);
  CHECK_STATUS_RET(*status, file);
  
  // Re-open the file.
  file=openImpactFile(filename, READWRITE, status);
  CHECK_STATUS_RET(*status, file);

  return(file);
}


void getNextImpactFromFile(ImpactFile* const file, Impact* const impact, 
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

  impact->energy = 0.;
  fits_read_col(file->fptr, TFLOAT, file->cenergy, file->row, 1, 1, 
		&impact->energy, &impact->energy, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->position.x = 0.;
  fits_read_col(file->fptr, TDOUBLE, file->cx, file->row, 1, 1, 
		&impact->position.x, &impact->position.x, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->position.y = 0.;
  fits_read_col(file->fptr, TDOUBLE, file->cy, file->row, 1, 1, 
		&impact->position.y, &impact->position.y, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->ph_id = 0;
  fits_read_col(file->fptr, TLONG, file->cph_id, file->row, 1, 1, 
		&impact->ph_id, &impact->ph_id, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->src_id = 0;
  fits_read_col(file->fptr, TLONG, file->csrc_id, file->row, 1, 1, 
		&impact->src_id, &impact->src_id, &anynul, status);
  CHECK_STATUS_VOID(*status);
  
  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("failed reading from impact list file");
    return;
  }

  return;
}


void addImpact2File(ImpactFile* const ilf, 
		    Impact* const impact, 
		    int* const status)
{
  ilf->row++;
  ilf->nrows++;

  fits_write_col(ilf->fptr, TDOUBLE, ilf->ctime, 
		 ilf->row, 1, 1, &impact->time, status);
  fits_write_col(ilf->fptr, TFLOAT, ilf->cenergy, 
		 ilf->row, 1, 1, &impact->energy, status);
  fits_write_col(ilf->fptr, TDOUBLE, ilf->cx, 
		 ilf->row, 1, 1, &(impact->position.x), status);
  fits_write_col(ilf->fptr, TDOUBLE, ilf->cy, 
		 ilf->row, 1, 1, &(impact->position.y), status);
  fits_write_col(ilf->fptr, TLONG, ilf->cph_id, 
		 ilf->row, 1, 1, &impact->ph_id, status);
  fits_write_col(ilf->fptr, TLONG, ilf->csrc_id, 
		 ilf->row, 1, 1, &impact->src_id, status);
}

