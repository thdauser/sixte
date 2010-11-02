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
  file->geneventfile=NULL;

  // Initialize.
  file->cpat_type=0;
  file->cphas    =0;
  file->geneventfile=newGenEventFile(status);

  return(file);
}



void destroyGenPatternFile(GenPatternFile** file, 
			   int* const status)
{
  if (NULL!=*file) {
    if (NULL!=(*file)->geneventfile) {
      destroyGenEventFile(&(*file)->geneventfile, status);
    }
    free(*file);
    *file=NULL;
  }
}



GenPatternFile* openNewGenPatternFile(const char* const filename,
				      const char* const template,
				      int* const status)
{
  GenPatternFile* file = newGenPatternFile(status);
  if (EXIT_SUCCESS!=*status) return(file);

  // Remove old file if it exists.
  remove(filename);

  // Open the GenEventFile.
  file->geneventfile = openNewGenEventFile(filename, template, status);

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
  file->geneventfile = openGenEventFile(filename, mode, status);

  if(fits_get_colnum(file->geneventfile->fptr, CASEINSEN, 
		     "PAT_TYPE", &file->cpat_type, status)) 
    return(file);
  if(fits_get_colnum(file->geneventfile->fptr, CASEINSEN, 
		     "PHAS", &file->cphas, status)) 
    return(file);

  return(file);
}



void addGenPattern2File(GenPatternFile* const file, 
			GenPattern* const pattern, 
			int* const status)
{
  // Call underlying routine.
  addGenEvent2File(file->geneventfile, &pattern->event, status);
  
  if (fits_write_col(file->geneventfile->fptr, TINT, 
		     file->cpat_type, file->geneventfile->row, 
		     1, 1, &pattern->pat_type, status)) return;
  if (fits_write_col(file->geneventfile->fptr, TLONG, 
		     file->cphas, file->geneventfile->row, 
		     1, 9, &pattern->phas, status)) return;
}



