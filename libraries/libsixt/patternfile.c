#include "patternfile.h"


PatternFile* newPatternFile(int* const status)
{
  PatternFile* file = (PatternFile*)malloc(sizeof(PatternFile));
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
  PatternFile* file = newPatternFile(status);
  if (EXIT_SUCCESS!=*status) return(file);

  // Remove old file if it exists.
  remove(filename);

  // Open the EventListFile.
  file->eventlistfile = openNewEventListFile(filename, status);

  // Close the file.
  destroyPatternFile(&file, status);
  if (EXIT_SUCCESS!=*status) return(file);

  // Re-open the file.
  file = openPatternFile(filename, READWRITE, status);
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

