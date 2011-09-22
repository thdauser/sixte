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
  file->cpat_type=0;
  file->cphas    =0;
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
		  "PAT_TYPE", &file->cpat_type, status);
  fits_get_colnum(file->eventlistfile->fptr, CASEINSEN, 
		  "PILEUP", &file->cpileup, status);
  fits_get_colnum(file->eventlistfile->fptr, CASEINSEN, 
		  "PHAS", &file->cphas, status);
  CHECK_STATUS_RET(*status, file);

  return(file);
}


void addPattern2File(PatternFile* const file, 
		     Pattern* const pattern, 
		     int* const status)
{
  // Call underlying routine.
  addEvent2File(file->eventlistfile, &pattern->event, status);
  
  fits_write_col(file->eventlistfile->fptr, TINT, 
		 file->cpat_type, file->eventlistfile->row, 
		 1, 1, &pattern->pat_type, status);
  fits_write_col(file->eventlistfile->fptr, TINT, 
		 file->cpileup, file->eventlistfile->row, 
		 1, 1, &pattern->pileup, status);
  fits_write_col(file->eventlistfile->fptr, TLONG, 
		 file->cphas, file->eventlistfile->row, 
		 1, 9, &pattern->phas, status);
  CHECK_STATUS_VOID(*status);
}

