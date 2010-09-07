#include "geneventfile.h"


GenEventFile* newGenEventFile(int* const status)
{
  GenEventFile* file = (GenEventFile*)malloc(sizeof(GenEventFile));
  if (NULL==file) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenEventFile failed!\n", 
		   *status);
    return(file);
  }

  // Initialize pointers with NULL.
  file->fptr=NULL;

  // Initialize values.
  file->nrows=0;
  file->row  =0;
  file->ctime=0;
  file->cpha =0;
  file->crawx=0;
  file->crawy=0;

  return(file);
}



void destroyGenEventFile(GenEventFile** file, int* const status)
{
  if (NULL!=*file) {
    if (NULL!=(*file)->fptr) {
      fits_close_file((*file)->fptr, status);
      (*file)->fptr=NULL;
    }
    free(*file);
    *file=NULL;
  }
}



GenEventFile* openNewGenEventFile(const char* const filename,
				  const char* const template,
				  int* const status)
{
  GenEventFile* file = newGenEventFile(status);
  if (EXIT_SUCCESS!=*status) return(file);

  // Remove old file if it exists.
  remove(filename);

  // Create a new event list FITS file from the template file.
  char buffer[MAXMSG];
  sprintf(buffer, "%s(%s)", filename, template);
  headas_chat(4, "create new event list file '%s' from template '%s' ...\n", 
	      filename, template);
  if (fits_create_file(&file->fptr, buffer, status)) return(file);

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
	if (fits_update_key(file->fptr, TSTRING, "DATE-OBS", current_time_str, 
			    "Start Time (UTC) of exposure", status)) 
	  return(file);
      }
    }
  } 
  // END of writing time information to Event File FITS header.

  // Add header information about program parameters.
  // The second parameter "1" means that the headers are written
  // to the first extension.
  HDpar_stamp(file->fptr, 1, status);

  // Move to the right (second) HDU with the binary table extension.
  int hdutype;
  if (fits_movabs_hdu(file->fptr, 2, &hdutype, status)) return(file);

  // Determine the column numbers.
  if(fits_get_colnum(file->fptr, CASEINSEN, "TIME", &file->ctime, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "PHA", &file->cpha, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "RAWX", &file->crawx, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "RAWY", &file->crawy, status)) 
    return(file);


  return(file);
}



void addGenEvent2File(GenEventFile* const file, GenEvent* const event, 
		      int* const status)
{
  // Check if the event file has been opened.
  if (NULL==file) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: no event file opened!\n", *status);
    return;
  }

  // Insert a new, empty row to the table:
  if (fits_insert_rows(file->fptr, file->row, 1, status)) return;
  file->row++;
  file->nrows++;

  // Insert the event data.
  if (fits_write_col(file->fptr, TDOUBLE, file->ctime, file->row, 
		     1, 1, &event->time, status)) return;
  if (fits_write_col(file->fptr, TLONG, file->cpha, file->row, 
		     1, 1, &event->pha, status)) return;
  if (fits_write_col(file->fptr, TINT, file->crawx, file->row, 
		     1, 1, &event->rawx, status)) return;
  if (fits_write_col(file->fptr, TINT, file->crawy, file->row, 
		     1, 1, &event->rawy, status)) return;

}

