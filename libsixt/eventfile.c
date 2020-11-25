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

#include "eventfile.h"


EventFile* newEventFile(int* const status)
{
  EventFile* file=(EventFile*)malloc(sizeof(EventFile));
  CHECK_NULL_RET(file, *status, "memory allocation for EventFile failed",
		 file);

  // Initialize pointers with NULL.
  file->fptr=NULL;
  file->memptr=NULL;

  // Initialize.
  file->nrows=0;
  file->memsize=0;
  file->ctime=0;
  file->cframe =0;
  file->cpha    =0;
  file->cpi    =0;
  file->csignal=0;
  file->crawx  =0;
  file->crawy  =0;
  file->cra    =0;
  file->cdec   =0;
  file->cph_id =0;
  file->csrc_id=0;
  file->ctype   =0;
  file->cnpixels=0;
  file->csignals=0;
  file->cphas    =0;
  file->cpileup =0;

  return(file);
}


void freeEventFile(EventFile** const file,
		   int* const status)
{
  if (NULL!=*file) {
    if (NULL!=(*file)->fptr) {
      // If the file was opened in READWRITE mode, calculate
      // the check sum an append it to the FITS header.
      int mode;
      fits_file_mode((*file)->fptr, &mode, status);
      if (READWRITE==mode) {
        fits_write_chksum((*file)->fptr, status);
      }
      fits_close_file((*file)->fptr, status);
    }

    free((*file)->memptr);
    free(*file);
    *file=NULL;
  }
}

void getColNumsFromEventFile(EventFile* const file,	int* const status)
{
	fits_get_colnum(file->fptr, CASEINSEN, "TIME", &file->ctime, status);
	fits_get_colnum(file->fptr, CASEINSEN, "FRAME", &file->cframe, status);
	fits_get_colnum(file->fptr, CASEINSEN, "PHA", &file->cpha, status);
	CHECK_STATUS_VOID(*status);
	fits_get_colnum(file->fptr, CASEINSEN, "PI", &file->cpi, status);
	if( *status == COL_NOT_FOUND ){   // Optional Column: reset status if not existent
		fits_clear_errmsg();
		*status=EXIT_SUCCESS;
	}
	fits_get_colnum(file->fptr, CASEINSEN, "SIGNAL", &file->csignal, status);
	fits_get_colnum(file->fptr, CASEINSEN, "RAWX", &file->crawx, status);
	fits_get_colnum(file->fptr, CASEINSEN, "RAWY", &file->crawy, status);
	fits_get_colnum(file->fptr, CASEINSEN, "RA", &file->cra, status);
	fits_get_colnum(file->fptr, CASEINSEN, "DEC", &file->cdec, status);
	fits_get_colnum(file->fptr, CASEINSEN, "PH_ID", &file->cph_id, status);
	fits_get_colnum(file->fptr, CASEINSEN, "SRC_ID", &file->csrc_id, status);
	fits_get_colnum(file->fptr, CASEINSEN, "NPIXELS", &file->cnpixels, status);
	fits_get_colnum(file->fptr, CASEINSEN, "TYPE", &file->ctype, status);
	fits_get_colnum(file->fptr, CASEINSEN, "PILEUP", &file->cpileup, status);
	fits_get_colnum(file->fptr, CASEINSEN, "SIGNALS", &file->csignals, status);
	fits_get_colnum(file->fptr, CASEINSEN, "PHAS", &file->cphas, status);
	CHECK_STATUS_VOID(*status);
}


EventFile* openNewEventFile(const char* const filename,
			    char* const telescop,
			    char* const instrume,
			    char* const filter,
			    char* const ancrfile,
			    char* const respfile,
			    const double mjdref,
			    const double timezero,
			    const double tstart,
			    const double tstop,
			    const int nxdim,
			    const int nydim,
			    const char clobber,
			    int* const status)
{
  fitsfile* fptr=NULL;
  CHECK_STATUS_RET(*status, NULL);

  // Check if the file already exists.
  int exists;
  fits_file_exists(filename, &exists, status);
  CHECK_STATUS_RET(*status, NULL);
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
      return(NULL);
    }
  }

  // Create a new event list FITS file from the template file.
  char buffer[MAXFILENAME];
  sprintf(buffer, "%s(%s%s)", filename, SIXT_DATA_PATH, "/templates/eventfile.tpl");
  fits_create_file(&fptr, buffer, status);
  if (*status!=EXIT_SUCCESS) {
      fits_report_error(stderr, *status);
      fprintf(stderr,"failed opening file %s\n",buffer);
  }
  CHECK_STATUS_RET(*status, NULL);

  // Insert header keywords to 1st HDU
  sixt_add_fits_stdkeywords_obsolete(fptr, 1, telescop, instrume, filter,
			    ancrfile, respfile,
			    mjdref, timezero, tstart, tstop, status);
  CHECK_STATUS_RET(*status, NULL);

  fits_update_key(fptr, TSTRING, "HDUCLAS1", "EVENTS", "", status);
  fits_update_key(fptr, TSTRING, "HDUCLASS", "OGIP", "", status);
  CHECK_STATUS_RET(*status, NULL);

  // Insert header keywords to 2nd HDU.
  sixt_add_fits_stdkeywords_obsolete(fptr, 2, telescop, instrume, filter,
			    ancrfile, respfile,
			    mjdref, timezero, tstart, tstop, status);
  CHECK_STATUS_RET(*status, NULL);

  // Close the file.
  fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, NULL);

  // Re-open the file.
  EventFile* plf=openEventFile(filename, READWRITE, status);
  CHECK_STATUS_RET(*status, plf);

  // Update the TLMIN and TLMAX keywords for the DETX and DETY columns.
  char keystr[MAXMSG];
  int ibuffer;
  sprintf(keystr, "TLMIN%d", plf->crawx);
  ibuffer=0;
  fits_update_key(plf->fptr, TINT, keystr, &ibuffer, "", status);
  sprintf(keystr, "TLMAX%d", plf->crawx);
  ibuffer=nxdim-1;
  fits_update_key(plf->fptr, TINT, keystr, &ibuffer, "", status);
  sprintf(keystr, "TLMIN%d", plf->crawy);
  ibuffer=0;
  fits_update_key(plf->fptr, TINT, keystr, &ibuffer, "", status);
  sprintf(keystr, "TLMAX%d", plf->crawy);
  ibuffer=nydim-1;
  fits_update_key(plf->fptr, TINT, keystr, &ibuffer, "", status);
  CHECK_STATUS_RET(*status, plf);

  return(plf);
}

EventFile* openEventFile(const char* const filename,
			 const int mode, int* const status)
{
  EventFile* file=newEventFile(status);
  CHECK_STATUS_RET(*status, file);

  headas_chat(3, "open event file '%s' ...\n", filename);
  fits_open_table(&file->fptr, filename, mode, status);
  CHECK_STATUS_RET(*status, file);

  // Determine the row numbers.
  fits_get_num_rows(file->fptr, &file->nrows, status);

  // Determine the column numbers.
  getColNumsFromEventFile(file, status);
  CHECK_STATUS_RET(*status, file);

  // Check if the vector length of the PH_ID and SRC_ID columns is equivalent
  // with the corresponding array lengths in the Event data structure.
  int typecode;
  long repeat, width;
  // PH_ID.
  fits_get_coltype(file->fptr, file->cph_id, &typecode, &repeat,
		   &width, status);
  CHECK_STATUS_RET(*status, file);
  if (repeat!=NEVENTPHOTONS) {
    // Throw an error.
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "inconsistent maximum number of photons contributing "
	    "to a single event (simulation: %d, event file template %ld)",
	    NEVENTPHOTONS, repeat);
    SIXT_ERROR(msg);
    return(file);
  }
  // SRC_ID.
  fits_get_coltype(file->fptr, file->csrc_id, &typecode, &repeat,
		   &width, status);
  CHECK_STATUS_RET(*status, file);
  if (repeat!=NEVENTPHOTONS) {
    // Throw an error.
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "inconsistent maximum number of photons contributing "
	    "to a single event (simulation: %d, event file template %ld)",
	    NEVENTPHOTONS, repeat);
    SIXT_ERROR(msg);
    return(file);
  }

  return(file);
}

void addCol2EventFile(EventFile* const file,
		int* const colnum,
		char* const ttype,
		char* const tform,
		char* const tunit,
		int* const status){

	// Check if colnum is out of bounds
	int cnum = 0;
	fits_get_num_cols(file->fptr, &cnum, status);
	CHECK_STATUS_VOID(*status);
	if( (*colnum <= 0) || (*colnum > cnum+1) ){
		*colnum = cnum+1;
	}

	// Insert column in File
	fits_insert_col(file->fptr,*colnum,ttype,tform,status);
	CHECK_STATUS_VOID(*status);

	// Set TUNIT of added column
	char keystr[MAXMSG];
	sprintf(keystr, "TUNIT%d", *colnum);
	char comment[MAXMSG];
	sprintf(comment, "Unit of column %s", ttype);
	fits_update_key(file->fptr, TSTRING, keystr, tunit, comment, status);
	CHECK_STATUS_VOID(*status);

	// Re-initilize column numbers in EventFile
	getColNumsFromEventFile(file, status);
	CHECK_STATUS_VOID(*status);

	// TEST if added column is really there
	fits_get_colnum(file->fptr, CASEINSEN, ttype, &cnum, status);
	if( *status == COL_NOT_FOUND ){
		char msg[MAXMSG];
    	sprintf(msg, "Adding column '%s' at %d failed!",
    			ttype,*colnum);
    	SIXT_ERROR(msg);
    	*status=EXIT_FAILURE;
	}
	return;
}

void addEvent2File(EventFile* const file,
		   Event* const event,
		   int* const status)
{
  // Check if the file has been opened.
  CHECK_NULL_VOID(file, *status, "event file not open");
  CHECK_NULL_VOID(file->fptr, *status, "event file not open");

  // Write the data.
  updateEventInFile(file, ++file->nrows, event, status);
  CHECK_STATUS_VOID(*status);
}


void getEventFromFile(const EventFile* const file,
		      const int row, Event* const event,
		      int* const status)
{
  // Check if the file has been opened.
  CHECK_NULL_VOID(file, *status, "event file not open");
  CHECK_NULL_VOID(file->fptr, *status, "event file not open");

  // Check if there is still a row available.
  if (row>file->nrows) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("event file contains no further entries");
    return;
  }

  // Read in the data.
  int anynul=0;
  double dnull=0.;
  float fnull=0.;
  long lnull=0;
  int inull=0;

  fits_read_col(file->fptr, TDOUBLE, file->ctime, row, 1, 1,
		&dnull, &event->time, &anynul, status);
  fits_read_col(file->fptr, TLONG, file->cframe, row, 1, 1,
		&lnull, &event->frame, &anynul, status);
  fits_read_col(file->fptr, TLONG, file->cpha, row, 1, 1,
		&lnull, &event->pha, &anynul, status);
  fits_read_col(file->fptr, TFLOAT, file->csignal, row, 1, 1,
		&fnull, &event->signal, &anynul, status);
  fits_read_col(file->fptr, TINT, file->crawx, row, 1, 1,
		&inull, &event->rawx, &anynul, status);
  fits_read_col(file->fptr, TINT, file->crawy, row, 1, 1,
		&inull, &event->rawy, &anynul, status);
  fits_read_col(file->fptr, TDOUBLE, file->cra, row, 1, 1,
		&dnull, &event->ra, &anynul, status);
  event->ra*=M_PI/180.;
  fits_read_col(file->fptr, TDOUBLE, file->cdec, row, 1, 1,
		&lnull, &event->dec, &anynul, status);
  event->dec*=M_PI/180.;
  fits_read_col(file->fptr, TLONG, file->cph_id, row, 1, NEVENTPHOTONS,
		&lnull, &event->ph_id, &anynul, status);
  fits_read_col(file->fptr, TLONG, file->csrc_id, row, 1, NEVENTPHOTONS,
		&lnull, &event->src_id, &anynul, status);
  fits_read_col(file->fptr, TLONG, file->cnpixels, row, 1, 1,
		&lnull, &event->npixels, &anynul, status);
  fits_read_col(file->fptr, TINT, file->ctype, row, 1, 1,
		&inull, &event->type, &anynul, status);
  fits_read_col(file->fptr, TINT, file->cpileup, row, 1, 1,
		&inull, &event->pileup, &anynul, status);
  fits_read_col(file->fptr, TFLOAT, file->csignals, row, 1, 9,
		&fnull, &event->signals, &anynul, status);
  fits_read_col(file->fptr, TLONG, file->cphas, row, 1, 9,
		&lnull, &event->phas, &anynul, status);
  CHECK_STATUS_VOID(*status);

  // only read PI column if file->cpi is valid
  if( file->cpi > 0 ){
  	fits_read_col(file->fptr, TLONG, file->cpi, row, 1, 1,
  			&dnull, &event->pi,&anynul, status);
  	CHECK_STATUS_VOID(*status);
  }

  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("reading from EventFile failed");
    return;
  }
}


void updateEventInFile(const EventFile* const file,
		       const int row, Event* const event,
		       int* const status)
{
//puts("write event.");
  fits_write_col(file->fptr, TDOUBLE, file->ctime, row,
		 1, 1, &event->time, status);
  fits_write_col(file->fptr, TLONG, file->cframe, row,
		 1, 1, &event->frame, status);
  fits_write_col(file->fptr, TLONG, file->cpha, row,
		 1, 1, &event->pha, status);
  fits_write_col(file->fptr, TFLOAT, file->csignal, row,
		 1, 1, &event->signal, status);
  fits_write_col(file->fptr, TINT, file->crawx, row,
		 1, 1, &event->rawx, status);
  fits_write_col(file->fptr, TINT, file->crawy, row,
		 1, 1, &event->rawy, status);
  double dbuffer=event->ra*180./M_PI;
  fits_write_col(file->fptr, TDOUBLE, file->cra, row,
		 1, 1, &dbuffer, status);
  dbuffer=event->dec*180./M_PI;
  fits_write_col(file->fptr, TDOUBLE, file->cdec, row,
		 1, 1, &dbuffer, status);
  fits_write_col(file->fptr, TLONG, file->cph_id, row,
		 1, NEVENTPHOTONS, &event->ph_id, status);
  fits_write_col(file->fptr, TLONG, file->csrc_id, row,
		 1, NEVENTPHOTONS, &event->src_id, status);
  fits_write_col(file->fptr, TLONG, file->cnpixels, row,
		 1, 1, &event->npixels, status);
  fits_write_col(file->fptr, TINT, file->cpileup, row,
		 1, 1, &event->pileup, status);
  fits_write_col(file->fptr, TINT, file->ctype, row,
		 1, 1, &event->type, status);
  fits_write_col(file->fptr, TFLOAT, file->csignals, row,
		 1, 9, &event->signals, status);
  fits_write_col(file->fptr, TLONG, file->cphas, row,
		 1, 9, &event->phas, status);
  CHECK_STATUS_VOID(*status);

  // only write PI value if event->pi value is valid, i.e., != -1.
  if( file->cpi > 0 ){
  	fits_write_col(file->fptr, TLONG, file->cpi, row,
  			1, 1, &event->pi, status);
    CHECK_STATUS_VOID(*status);
  }
}


void copyEventFile(const EventFile* const src,
		   EventFile* const dest,
		   const float threshold_lo_keV,
		   const float threshold_up_keV,
		   int* const status)
{
  // Check if the event file is empty.
  if (dest->nrows>0) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("destination event file is not empty");
    return;
  }

  // Copy the event type.
  char evtype[MAXMSG], comment[MAXMSG];
  fits_read_key(src->fptr, TSTRING, "EVTYPE", evtype, comment, status);
  if (EXIT_SUCCESS!=*status) {
    SIXT_ERROR("could not read FITS keyword 'EVTYPE'");
    return;
  }
  fits_update_key(dest->fptr, TSTRING, "EVTYPE", evtype, comment, status);
  CHECK_STATUS_VOID(*status);

  // Get memory for buffers.
  Event* event=getEvent(status);
  CHECK_STATUS_VOID(*status);

  // Loop over all rows in the event file.
  long row;
  for (row=0; row<src->nrows; row++) {

    // Read an event from the input list.
    getEventFromFile(src, row+1, event, status);
    CHECK_STATUS_BREAK(*status);

    // Apply the lower event threshold.
    if (event->signal<threshold_lo_keV) {
      continue;
    }

    // Apply the upper event threshold.
    if ((threshold_up_keV>0.0)&&(event->signal>threshold_up_keV)) {
      continue;
    }

    // Add the new event to the output file.
    addEvent2File(dest, event, status);
    CHECK_STATUS_BREAK(*status);
  }
  CHECK_STATUS_VOID(*status);

  // Free memory.
  freeEvent(&event);
}

EventFile* copyEventFileMemory(EventFile* infile, int* const status) {
  // Create a new eventfile (fptr not initialized yet)
  EventFile* outfile = newEventFile(status);

  // Copy housekeeping data
  outfile->nrows = infile->nrows;
  outfile->ctime = infile->ctime;
  outfile->cframe = infile->cframe;
  outfile->cpha = infile->cpha;
  outfile->cpi = infile->cpi;
  outfile->csignal = infile->csignal;
  outfile->crawx = infile->crawx;
  outfile->crawy = infile->crawy;
  outfile->cra = infile->cra;
  outfile->cdec = infile->cdec;
  outfile->cph_id = infile->cph_id;
  outfile->csrc_id = infile->csrc_id;
  outfile->cnpixels = infile->cnpixels;
  outfile->ctype = infile->ctype;
  outfile->cpileup = infile->cpileup;
  outfile->csignals = infile->csignals;
  outfile->cphas = infile->cphas;

  // Create new fitsfile in core memory
  size_t deltasize = 1000000; // TODO: hardcoded for now; perhaps these
  outfile->memsize = 2880;    // values can be optimized.
  outfile->memptr = malloc(outfile->memsize);
  fits_create_memfile(&(outfile->fptr), &(outfile->memptr), &(outfile->memsize), deltasize,
                      realloc, status);

  // Copy all HDUs from infile to outfile
  fits_copy_file(infile->fptr, outfile->fptr, 1, 1, 1, status);

  return outfile;
}
