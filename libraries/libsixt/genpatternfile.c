#include "genpatternfile.h"


GenPatternFile* newGenPatternFile(int* const status)
{
  GenPatternFile* file = (GenPatternFile*)malloc(sizeof(GenPatternFile));
  if (NULL==file) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenPatternFile failed!\n", 
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



void destroyGenPatternFile(GenPatternFile** const file, 
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



GenPatternFile* openNewGenPatternFile(const char* const filename,
				      const char* const template,
				      double mjdref,
				      int* const status)
{
  GenPatternFile* file = newGenPatternFile(status);
  if (EXIT_SUCCESS!=*status) return(file);

  // Remove old file if it exists.
  remove(filename);

  // Open the EventListFile.
  file->eventlistfile = 
    openNewEventListFile(filename, template, mjdref, status);

  // Close the file.
  destroyGenPatternFile(&file, status);
  if (EXIT_SUCCESS!=*status) return(file);

  // Re-open the file.
  file = openGenPatternFile(filename, READWRITE, status);
  if (EXIT_SUCCESS!=*status) return(file);
  
  return(file);
}



GenPatternFile* openGenPatternFile(const char* const filename,
				   const int mode, int* const status)
{
  GenPatternFile* file = newGenPatternFile(status);
  if (EXIT_SUCCESS!=*status) return(file);

  headas_chat(4, "open pattern file '%s' ...\n", filename);

  // Call underlying open routine.
  file->eventlistfile = openEventListFile(filename, mode, status);

  if(fits_get_colnum(file->eventlistfile->fptr, CASEINSEN, 
		     "PAT_TYPE", &file->cpat_type, status)) 
    return(file);
  if(fits_get_colnum(file->eventlistfile->fptr, CASEINSEN, 
		     "PILEUP", &file->cpileup, status)) 
    return(file);
  if(fits_get_colnum(file->eventlistfile->fptr, CASEINSEN, 
		     "PHAS", &file->cphas, status)) 
    return(file);

  return(file);
}



void addGenPattern2File(GenPatternFile* const file, 
			GenPattern* const pattern, 
			int* const status)
{
  // Call underlying routine.
  addEvent2File(file->eventlistfile, &pattern->event, status);
  
  if (fits_write_col(file->eventlistfile->fptr, TINT, 
		     file->cpat_type, file->eventlistfile->row, 
		     1, 1, &pattern->pat_type, status)) return;
  if (fits_write_col(file->eventlistfile->fptr, TINT, 
		     file->cpileup, file->eventlistfile->row, 
		     1, 1, &pattern->pileup, status)) return;
  if (fits_write_col(file->eventlistfile->fptr, TLONG, 
		     file->cphas, file->eventlistfile->row, 
		     1, 9, &pattern->phas, status)) return;
}



