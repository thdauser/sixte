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


   Copyright 2014 Thorsten Brand, FAU
*/

#include "pixelimpactfile.h"

PixImpFile* newPixImpFile(int* const status){
  
  PixImpFile* file=(PixImpFile*)malloc(sizeof(PixImpFile));
  if (NULL==file) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for PixImpFile failed");
    return(file);
  }

  // Initialize pointers with NULL.
  file->fptr=NULL;

  // Initialize values.
  file->nrows	=0;
  file->row  	=0;
  file->ctime	=0;
  file->cenergy	=0;
  file->cph_id 	=0;
  file->csrc_id	=0;
  file->cpix_id	=0;
  file->cu	=0;
  file->cv	=0;
  file->cx	=0;
  file->cy	=0;
  file->cgrade1 = 0;
  file->cgrade2 = 0;
  file->ctotalenergy = 0;

  return(file);
}

void freePixImpFile(PixImpFile** const file, int* const status){
  
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

PixImpFile* openPixImpFile(const char* const filename,
			   const int mode, int* const status)
{
  PixImpFile* file=newPixImpFile(status);
  CHECK_STATUS_RET(*status, NULL);
  
  headas_chat(5, "open pixel impact list file '%s' ...\n", filename);
  
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
    return(NULL);
  }
  
  // Determine the number of rows in the impact list.
  if (fits_get_num_rows(file->fptr, &file->nrows, status)) return(file);
  // Set internal row counter to the beginning of the file (starting at 0).
  file->row = 0;

  // Determine the individual column numbers.
  fits_get_colnum(file->fptr, CASEINSEN, "TIME", &file->ctime, status);
  fits_get_colnum(file->fptr, CASEINSEN, "ENERGY", &file->cenergy, status);
  fits_get_colnum(file->fptr, CASEINSEN, "U", &file->cu, status);
  fits_get_colnum(file->fptr, CASEINSEN, "V", &file->cv, status);
  fits_get_colnum(file->fptr, CASEINSEN, "X", &file->cx, status);
  fits_get_colnum(file->fptr, CASEINSEN, "Y", &file->cy, status);
  fits_get_colnum(file->fptr, CASEINSEN, "PH_ID", &file->cph_id, status);
  fits_get_colnum(file->fptr, CASEINSEN, "SRC_ID", &file->csrc_id, status);
  fits_get_colnum(file->fptr, CASEINSEN, "PIXID", &file->cpix_id, status);
  CHECK_STATUS_RET(*status, NULL);
  fits_get_colnum(file->fptr, CASEINSEN, "GRADE1", &file->cgrade1, status);
  if (*status==COL_NOT_FOUND) {
    file->cgrade1=-1;
    *status=0;
  }
  CHECK_STATUS_RET(*status, NULL);
  fits_get_colnum(file->fptr, CASEINSEN, "GRADE2", &file->cgrade2, status);
  if (*status==COL_NOT_FOUND) {
    file->cgrade2=-1;
    *status=0;
  }
  CHECK_STATUS_RET(*status, NULL);
  fits_get_colnum(file->fptr, CASEINSEN, "TOTALEN", &file->ctotalenergy, status);
  if (*status==COL_NOT_FOUND) {
    file->ctotalenergy=-1;
    *status=0;
  }
  CHECK_STATUS_RET(*status, NULL);
  return(file);
}

PixImpFile* openNewPixImpFile(const char* const filename,
			      char* const telescop,
			      char* const instrume,
			      char* const filter,
			      char* const ancrfile,
			      char* const respfile,
			      char* const xmlfile,
			      char* const impactlist,
			      const double mjdref,
			      const double timezero,
			      const double tstart,
			      const double tstop,
			      const char clobber,
			      int* const status){
  
  PixImpFile* file=newPixImpFile(status);
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
	  "/templates/pixelimpactfile_gradingcols.tpl");
  fits_create_file(&file->fptr, buffer, status);
  CHECK_STATUS_RET(*status, file);

  // Insert header keywords to 1st and 2nd HDU.
  sixt_add_fits_stdkeywords_obsolete(file->fptr, 1, telescop, instrume, filter,
			    ancrfile, respfile,
			    mjdref, timezero, tstart, tstop, status);
  CHECK_STATUS_RET(*status, NULL);
  sixt_add_fits_stdkeywords_obsolete(file->fptr, 2, telescop, instrume, filter,
			    ancrfile, respfile,
			    mjdref, timezero, tstart, tstop, status);
  CHECK_STATUS_RET(*status, NULL);
  
  // Move to the binary table extension.
  fits_movabs_hdu(file->fptr, 2, 0, status);
  CHECK_STATUS_RET(*status, file);
  
  //Write XML into header
  char comment[MAXMSG];
  sprintf(comment, "XMLFILE: %s", xmlfile);
  fits_write_comment(file->fptr, comment, status);
  CHECK_STATUS_RET(*status, file);
  
  //Write photon impact file into header
  sprintf(comment, "PHOFILE: %s", impactlist);
  fits_write_comment(file->fptr, comment, status);
  CHECK_STATUS_RET(*status, file);

  // Close the new ImpactFile.
  freePixImpFile(&file, status);
  CHECK_STATUS_RET(*status, file);
  
  // Re-open the file.
  file=openPixImpFile(filename, READWRITE, status);
  CHECK_STATUS_RET(*status, file);

  return(file);
}

int getNextImpactFromPixImpFile(PixImpFile* const file, 
			   PixImpact* const impact,
			   int* const status)
{
  // Check if the file has been opened.
  if (NULL==file) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("no pixel impact list file opened");
    return 0;
  }
  if (NULL==file->fptr) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("no pixel impact list file opened");
    return 0;
  }

  // Move counter to next line.
  file->row++;

  // Check if there is still a row available.
  if (file->row > file->nrows) {
    return 0;
  }

  // Read in the data.
  // JW: note: initialization of the data is NOT necessary 
  // (only costs time)
  int anynul=0;
  fits_read_col(file->fptr, TDOUBLE, file->ctime, file->row, 1, 1, 
		NULL, &impact->time, &anynul, status);
  fits_read_col(file->fptr, TFLOAT, file->cenergy, file->row, 1, 1, 
		NULL, &impact->energy, &anynul, status);
  fits_read_col(file->fptr, TDOUBLE, file->cx, file->row, 1, 1, 
		NULL, &impact->detposition.x, &anynul, status);
  fits_read_col(file->fptr, TDOUBLE, file->cy, file->row, 1, 1, 
		NULL, &impact->detposition.y, &anynul, status);
  fits_read_col(file->fptr, TDOUBLE, file->cu, file->row, 1, 1, 
		NULL, &impact->pixposition.x, &anynul, status);
  fits_read_col(file->fptr, TDOUBLE, file->cv, file->row, 1, 1, 
		NULL, &impact->pixposition.y, &anynul, status);
  fits_read_col(file->fptr, TLONG, file->cph_id, file->row, 1, 1, 
		NULL, &impact->ph_id, &anynul, status);
  fits_read_col(file->fptr, TLONG, file->csrc_id, file->row, 1, 1, 
		NULL, &impact->src_id, &anynul, status);
  fits_read_col(file->fptr, TLONG, file->cpix_id, file->row, 1, 1, 
		NULL, &impact->pixID, &anynul, status);
  impact->pixID=impact->pixID-1; // why??
  
  if (file->cgrade1!=-1) {
    fits_read_col(file->fptr, TLONG, file->cgrade1, file->row, 1, 1,
		  NULL, &impact->grade1, &anynul, status);
  }
  if (file->cgrade2!=-1) {
    fits_read_col(file->fptr, TLONG, file->cgrade2, file->row, 1, 1,
		  NULL, &impact->grade2, &anynul, status);
  }
  if (file->ctotalenergy!=-1) {
    fits_read_col(file->fptr, TDOUBLE, file->ctotalenergy, file->row, 1, 1,
		  NULL, &impact->totalenergy, &anynul, status);
  }
  
  // Check if an error occurred during the reading process.
  if (*status!=0) {
    SIXT_ERROR("failed reading from pixel impact list file");
    return 0;
  }

  return 1;
}

void addImpact2PixImpFile(PixImpFile* const ilf, 
			  PixImpact* const impact,
			  int* const status)
{
  ilf->row++;
  ilf->nrows++;
  
  impact->pixID=impact->pixID+1;

  fits_write_col(ilf->fptr, TDOUBLE, ilf->ctime, 
		 ilf->row, 1, 1, &impact->time, status);
  fits_write_col(ilf->fptr, TFLOAT, ilf->cenergy, 
		 ilf->row, 1, 1, &impact->energy, status);
  fits_write_col(ilf->fptr, TDOUBLE, ilf->cu, 
		 ilf->row, 1, 1, &(impact->pixposition.x), status);
  fits_write_col(ilf->fptr, TDOUBLE, ilf->cv, 
		 ilf->row, 1, 1, &(impact->pixposition.y), status);
  fits_write_col(ilf->fptr, TDOUBLE, ilf->cx, 
		 ilf->row, 1, 1, &(impact->detposition.x), status);
  fits_write_col(ilf->fptr, TDOUBLE, ilf->cy, 
		 ilf->row, 1, 1, &(impact->detposition.y), status);
  fits_write_col(ilf->fptr, TLONG, ilf->cph_id, 
		 ilf->row, 1, 1, &impact->ph_id, status);
  fits_write_col(ilf->fptr, TLONG, ilf->csrc_id, 
		 ilf->row, 1, 1, &impact->src_id, status);
  fits_write_col(ilf->fptr, TLONG, ilf->cpix_id, 
		 ilf->row, 1, 1, &impact->pixID, status);

  long default_grade=0;
  double default_totalenergy = 0.;
  if (ilf->cgrade1!=-1) {
	  fits_write_col(ilf->fptr, TLONG, ilf->cgrade1,
			 ilf->row, 1, 1, &default_grade, status);
  }
  if (ilf->cgrade2!=-1) {
	  fits_write_col(ilf->fptr, TLONG, ilf->cgrade2,
			 ilf->row, 1, 1, &default_grade, status);
  }
  if (ilf->ctotalenergy!=-1) {
	  fits_write_col(ilf->fptr, TDOUBLE, ilf->ctotalenergy,
			 ilf->row, 1, 1, &default_totalenergy, status);
  }
}

void updateGradingPixImp(PixImpFile* const ilf,
		  int row, long grade1, long grade2,double totalenergy,
		  int* const status){
	if (row>ilf->nrows){
		*status=EXIT_FAILURE;
		SIXT_ERROR("updating of grading in pixImpact File failed");
		return;

	}
	fits_write_col(ilf->fptr, TLONG, ilf->cgrade1,
			row, 1, 1, &grade1, status);
	fits_write_col(ilf->fptr, TLONG, ilf->cgrade2,
			row, 1, 1, &grade2, status);
	fits_write_col(ilf->fptr, TDOUBLE, ilf->ctotalenergy,
			row, 1, 1, &totalenergy, status);

}

void getPixImpFileTimeValues(PixImpFile* const ilf,
			  double *mjdref,
			  double *timezero,
			  double *tstart,
			  double *tstop,
			  int* const status)
{
  char comment[MAXMSG];
  fits_read_key(ilf->fptr, TDOUBLE, "MJDREF", mjdref, comment, status);
  CHECK_STATUS_VOID(*status);
  fits_read_key(ilf->fptr, TDOUBLE, "TIMEZERO", timezero, comment, status);
  CHECK_STATUS_VOID(*status);
  fits_read_key(ilf->fptr, TDOUBLE, "TSTART", tstart, comment, status);
  CHECK_STATUS_VOID(*status);
  fits_read_key(ilf->fptr, TDOUBLE, "TSTOP", tstop, comment, status);
  CHECK_STATUS_VOID(*status);
}

void copyPixImpact(PixImpact* impact1, PixImpact* impact2){
	impact1->pixID = impact2->pixID;
	impact1->time = impact2->time;
	impact1->energy = impact2->energy;
	impact1->detposition.x = impact2->detposition.x;
	impact1->detposition.y = impact2->detposition.y;
	impact1->pixposition.x = impact2->pixposition.x;
	impact1->pixposition.y = impact2->pixposition.y;
	impact1->ph_id = impact2->ph_id;
	impact1->src_id = impact2->src_id;
	impact1->nb_pileup = impact2->nb_pileup;
	impact1->weight_index=impact2->weight_index;
}
