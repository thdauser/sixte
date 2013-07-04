#include "patternfile.h"


PatternFile* newPatternFile(int* const status)
{
  PatternFile* file=(PatternFile*)malloc(sizeof(PatternFile));
  CHECK_NULL_RET(file, *status, "memory allocation for PatternFile failed", 
		 file);

  // Initialize pointers with NULL.
  file->fptr=NULL;

  // Initialize.
  file->nrows=0;
  file->ctime=0;
  file->cframe =0;
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
  file->cpis    =0;
  file->cpileup =0;

  return(file);
}


void destroyPatternFile(PatternFile** const file, 
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
    free(*file);
    *file=NULL;
  }
}


PatternFile* openNewPatternFile(const char* const filename,
				char* const telescop,
				char* const instrume,
				char* const filter,
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
  sprintf(buffer, "%s(%s%s)", filename, SIXT_DATA_PATH, 
	  "/templates/patternlist.tpl");
  fits_create_file(&fptr, buffer, status);
  CHECK_STATUS_RET(*status, NULL);

  // Insert header keywords to 1st and 2nd HDU.
  sixt_add_fits_stdkeywords(fptr, 1, telescop, instrume, filter,
			    mjdref, timezero, tstart, tstop, status);
  CHECK_STATUS_RET(*status, NULL);
  sixt_add_fits_stdkeywords(fptr, 2, telescop, instrume, filter,
			    mjdref, timezero, tstart, tstop, status);
  CHECK_STATUS_RET(*status, NULL);

  // Close the file.
  fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, NULL);

  // Re-open the file.
  PatternFile* plf=openPatternFile(filename, READWRITE, status);
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


PatternFile* openPatternFile(const char* const filename,
			     const int mode, int* const status)
{
  PatternFile* file = newPatternFile(status);
  CHECK_STATUS_RET(*status, file);

  headas_chat(4, "open pattern file '%s' ...\n", filename);
  fits_open_table(&file->fptr, filename, mode, status);
  CHECK_STATUS_RET(*status, file);

  // Determine the row numbers.
  fits_get_num_rows(file->fptr, &file->nrows, status);

  // Determine the column numbers.
  fits_get_colnum(file->fptr, CASEINSEN, "TIME", &file->ctime, status);
  fits_get_colnum(file->fptr, CASEINSEN, "FRAME", &file->cframe, status);
  fits_get_colnum(file->fptr, CASEINSEN, "PI", &file->cpi, status);
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
  fits_get_colnum(file->fptr, CASEINSEN, "PIS", &file->cpis, status);
  CHECK_STATUS_RET(*status, file);

  // Check if the vector length of the PH_ID and SRC_ID columns is equivalent 
  // with the corresponding array lengths in the Event data structure.
  int typecode;
  long repeat, width;
  // PH_ID.
  fits_get_coltype(file->fptr, file->cph_id, &typecode, &repeat,
		   &width, status);
  CHECK_STATUS_RET(*status, file);
  if (repeat!=NPATTERNPHOTONS) {
    // Throw an error.
    *status = EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "Error: the maximum number of photons contributing "
	    "to a single pattern is different for the parameter set "
	    "in the simulation (%d) and in the pattern list "
	    "template file (%ld)!\n", NPATTERNPHOTONS, repeat);
    HD_ERROR_THROW(msg, *status);
    return(file);
  }
  // SRC_ID.
  fits_get_coltype(file->fptr, file->csrc_id, &typecode, &repeat,
		   &width, status);
  CHECK_STATUS_RET(*status, file);
  if (repeat!=NPATTERNPHOTONS) {
    // Throw an error.
    *status = EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "Error: the maximum number of photons contributing "
	    "to a single pattern is different for the parameter set "
	    "in the simulation (%d) and in the pattern list "
	    "template file (%ld)!\n", NPATTERNPHOTONS, repeat);
    HD_ERROR_THROW(msg, *status);
    return(file);
  }

  return(file);
}


void addPattern2File(PatternFile* const file, 
		     Pattern* const pattern, 
		     int* const status)
{
  // Check if the file has been opened.
  CHECK_NULL_VOID(file, *status, "pattern file not open");
  CHECK_NULL_VOID(file->fptr, *status, "pattern file not open");

  // Write the data.
  updatePatternInFile(file, ++file->nrows, pattern, status);
  CHECK_STATUS_VOID(*status);
}


void getPatternFromFile(const PatternFile* const file,
			const int row, Pattern* const pattern,
			int* const status)
{
  // Check if the file has been opened.
  CHECK_NULL_VOID(file, *status, "pattern file not open");
  CHECK_NULL_VOID(file->fptr, *status, "pattern file not open");

  // Check if there is still a row available.
  if (row > file->nrows) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("pattern file contains no further entries");
    return;
  }

  // Read in the data.
  int anynul=0;
  double dnull=0.;
  float fnull=0.;
  long lnull=0;
  int inull=0;

  fits_read_col(file->fptr, TDOUBLE, file->ctime, row, 1, 1, 
		&dnull, &pattern->time, &anynul, status);
  fits_read_col(file->fptr, TLONG, file->cframe, row, 1, 1, 
		&lnull, &pattern->frame, &anynul, status);
  fits_read_col(file->fptr, TLONG, file->cpi, row, 1, 1, 
		&lnull, &pattern->pi, &anynul, status);
  fits_read_col(file->fptr, TFLOAT, file->csignal, row, 1, 1, 
		&fnull, &pattern->signal, &anynul, status);
  fits_read_col(file->fptr, TINT, file->crawx, row, 1, 1, 
		&inull, &pattern->rawx, &anynul, status);
  fits_read_col(file->fptr, TINT, file->crawy, row, 1, 1, 
		&inull, &pattern->rawy, &anynul, status);
  fits_read_col(file->fptr, TDOUBLE, file->cra, row, 1, 1, 
		&dnull, &pattern->ra, &anynul, status);
  pattern->ra *= M_PI/180.;
  fits_read_col(file->fptr, TDOUBLE, file->cdec, row, 1, 1, 
		&lnull, &pattern->dec, &anynul, status);
  pattern->dec*= M_PI/180.;
  fits_read_col(file->fptr, TLONG, file->cph_id, row, 1, NPATTERNPHOTONS, 
		&lnull, &pattern->ph_id, &anynul, status);
  fits_read_col(file->fptr, TLONG, file->csrc_id, row, 1, NPATTERNPHOTONS, 
		&lnull, &pattern->src_id, &anynul, status);
  fits_read_col(file->fptr, TLONG, file->cnpixels, row, 1, 1, 
		&lnull, &pattern->npixels, &anynul, status);
  fits_read_col(file->fptr, TINT, file->ctype, row, 1, 1, 
		&inull, &pattern->type, &anynul, status);
  fits_read_col(file->fptr, TINT, file->cpileup, row, 1, 1, 
		&inull, &pattern->pileup, &anynul, status);
  fits_read_col(file->fptr, TFLOAT, file->csignals, row, 1, 9, 
		&fnull, &pattern->signals, &anynul, status);
  fits_read_col(file->fptr, TLONG, file->cpis, row, 1, 9, 
		&lnull, &pattern->pis, &anynul, status);
  CHECK_STATUS_VOID(*status);

  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("reading from PatternFile failed");
    return;
  }
}


void updatePatternInFile(const PatternFile* const file,
			 const int row, Pattern* const pattern,
			 int* const status)
{
  fits_write_col(file->fptr, TDOUBLE, file->ctime, row, 
		 1, 1, &pattern->time, status);
  fits_write_col(file->fptr, TLONG, file->cframe, row, 
		 1, 1, &pattern->frame, status);
  fits_write_col(file->fptr, TLONG, file->cpi, row, 
		 1, 1, &pattern->pi, status);
  fits_write_col(file->fptr, TFLOAT, file->csignal, row, 
		 1, 1, &pattern->signal, status);
  fits_write_col(file->fptr, TINT, file->crawx, row, 
		 1, 1, &pattern->rawx, status);
  fits_write_col(file->fptr, TINT, file->crawy, row, 
		 1, 1, &pattern->rawy, status);
  double dbuffer = pattern->ra  * 180./M_PI;
  fits_write_col(file->fptr, TDOUBLE, file->cra, row, 
		 1, 1, &dbuffer, status);
  dbuffer = pattern->dec * 180./M_PI;
  fits_write_col(file->fptr, TDOUBLE, file->cdec, row, 
		 1, 1, &dbuffer, status);
  fits_write_col(file->fptr, TLONG, file->cph_id, row, 
		 1, NPATTERNPHOTONS, &pattern->ph_id, status);
  fits_write_col(file->fptr, TLONG, file->csrc_id, row, 
		 1, NPATTERNPHOTONS, &pattern->src_id, status);
  fits_write_col(file->fptr, TLONG, file->cnpixels, row, 
		 1, 1, &pattern->npixels, status);
  fits_write_col(file->fptr, TINT, file->cpileup, row, 
		 1, 1, &pattern->pileup, status);
  fits_write_col(file->fptr, TINT, file->ctype, row, 
		 1, 1, &pattern->type, status);
  fits_write_col(file->fptr, TFLOAT, file->csignals, row, 
		 1, 9, &pattern->signals, status);
  fits_write_col(file->fptr, TLONG, file->cpis, row, 
		 1, 9, &pattern->pis, status);
  CHECK_STATUS_VOID(*status);
}


void copyEvents2PatternFile(const EventListFile* const elf,
			    PatternFile* const plf,
			    int* const status)
{
  // Check if the pattern file is empty.
  if (plf->nrows>0) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("pattern file is not empty");
    return;
  }
  
  // Get memory for buffers.
  Event* event=getEvent(status);
  CHECK_STATUS_VOID(*status);
  Pattern* pattern=getPattern(status);
  CHECK_STATUS_VOID(*status);

  // Loop over all rows in the event file.
  long row;
  for (row=0; row<elf->nrows; row++) {

    // Read an event from the input list.
    getEventFromFile(elf, row+1, event, status);
    CHECK_STATUS_BREAK(*status);
    
    // Copy event data to pattern.
    pattern->rawx   =event->rawx;
    pattern->rawy   =event->rawy;
    pattern->time   =event->time;
    pattern->frame  =event->frame;
    pattern->pi     =event->pi;
    pattern->signal =event->signal;
    pattern->ra     =0.;
    pattern->dec    =0.;
    pattern->npixels=1;
    pattern->type   =0;
    
    pattern->pileup =0;
    int ii;
    for (ii=0; (ii<NEVENTPHOTONS)&&(ii<NPATTERNPHOTONS); ii++){
      pattern->ph_id[ii] =event->ph_id[ii];
      pattern->src_id[ii]=event->src_id[ii];

      if ((ii>0)&&(pattern->ph_id[ii]!=0)) {
	pattern->pileup=1;
      }
    }
    
    pattern->signals[0]=0.;
    pattern->signals[1]=0.;
    pattern->signals[2]=0.;
    pattern->signals[3]=0.;
    pattern->signals[4]=event->signal;
    pattern->signals[5]=0.;
    pattern->signals[6]=0.;
    pattern->signals[7]=0.;
    pattern->signals[8]=0.;

    // Add the new pattern to the output file.
    addPattern2File(plf, pattern, status);	  
    CHECK_STATUS_BREAK(*status);
  }
  CHECK_STATUS_VOID(*status);

  // Free memory.
  freeEvent(&event);
  freePattern(&pattern);
}

