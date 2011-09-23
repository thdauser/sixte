#include "patternfile.h"


PatternFile* newPatternFile(int* const status)
{
  PatternFile* file=(PatternFile*)malloc(sizeof(PatternFile));
  if (NULL==file) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for PatternFile failed!\n", 
		   *status);
    return(file);
  }

  // Initialize pointers with NULL.
  file->eventlistfile=NULL;

  // Initialize.
  file->ctype    =0;
  file->cnpixels =0;
  file->csignals =0;
  file->cpileup  =0;

  return(file);
}


void destroyPatternFile(PatternFile** const file, 
			int* const status)
{
  if (NULL!=*file) {
    if (NULL!=(*file)->eventlistfile) {
      freeEventListFile(&(*file)->eventlistfile, status);
    }
    free(*file);
    *file=NULL;
  }
}


PatternFile* openNewPatternFile(const char* const filename,
				int* const status)
{
  PatternFile* file=NULL;

  // Remove old file if it exists.
  remove(filename);

  // Create a new event list FITS file from the template file.
  fitsfile* fptr=NULL;
  char buffer[MAXFILENAME];
  sprintf(buffer, "%s(%s%s)", filename, SIXT_DATA_PATH, "/templates/patternlist.tpl");
  fits_create_file(&fptr, buffer, status);
  CHECK_STATUS_RET(*status, file);

  // Set the time-keyword in the Event List Header.
  // See also: Stevens, "Advanced Programming in the UNIX environment",
  // p. 155 ff.
  time_t current_time;
  if (0 != time(&current_time)) {
    struct tm* current_time_utc = gmtime(&current_time);
    if (NULL != current_time_utc) {
      char current_time_str[MAXMSG];
      if (strftime(current_time_str, MAXMSG, "%Y-%m-%dT%H:%M:%S", 
		   current_time_utc) > 0) {
	// Return value should be == 19 !
	if (fits_update_key(fptr, TSTRING, "DATE-OBS", current_time_str, 
			    "Start Time (UTC) of exposure", status)) 
	  return(file);
      }
    }
  } 
  // END of writing time information to pattern file FITS header.

  // Add header information about program parameters.
  // The second parameter "1" means that the headers are written
  // to the first extension.
  HDpar_stamp(fptr, 1, status);
  if (EXIT_SUCCESS!=*status) return(file);

  // Close the file.
  fits_close_file(fptr, status);
  if (EXIT_SUCCESS!=*status) return(file);

  // Re-open the file.
  file=openPatternFile(filename, READWRITE, status);
  if (EXIT_SUCCESS!=*status) return(file);
  
  return(file);
}


PatternFile* openPatternFile(const char* const filename,
			     const int mode, int* const status)
{
  PatternFile* file = newPatternFile(status);
  CHECK_STATUS_RET(*status, file);

  headas_chat(4, "open pattern file '%s' ...\n", filename);

  // Call underlying open routine.
  file->eventlistfile = openEventListFile(filename, mode, status);
  CHECK_STATUS_RET(*status, file);

  fits_get_colnum(file->eventlistfile->fptr, CASEINSEN, 
		  "NPIXELS", &file->cnpixels, status);
  fits_get_colnum(file->eventlistfile->fptr, CASEINSEN, 
		  "TYPE", &file->ctype, status);
  fits_get_colnum(file->eventlistfile->fptr, CASEINSEN, 
		  "PILEUP", &file->cpileup, status);
  fits_get_colnum(file->eventlistfile->fptr, CASEINSEN, 
		  "SIGNALS", &file->csignals, status);
  CHECK_STATUS_RET(*status, file);

  return(file);
}


void addPattern2File(PatternFile* const file, 
		     Pattern* const pattern, 
		     int* const status)
{
  // Call underlying routine.
  addEvent2File(file->eventlistfile, pattern->event, status);
  
  fits_write_col(file->eventlistfile->fptr, TLONG, 
		 file->cnpixels, file->eventlistfile->row, 
		 1, 1, &pattern->npixels, status);
  fits_write_col(file->eventlistfile->fptr, TINT, 
		 file->ctype, file->eventlistfile->row, 
		 1, 1, &pattern->type, status);
  fits_write_col(file->eventlistfile->fptr, TINT, 
		 file->cpileup, file->eventlistfile->row, 
		 1, 1, &pattern->pileup, status);
  fits_write_col(file->eventlistfile->fptr, TFLOAT, 
		 file->csignals, file->eventlistfile->row, 
		 1, 9, &pattern->signals, status);
  CHECK_STATUS_VOID(*status);
}

