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
  file->cframe =0;
  file->cpileup=0;
  file->cpat_type=0;
  file->cpat_id  =0;
  file->cpat_alig=0;

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
  if (EXIT_SUCCESS!=*status) return(file);

  // Close the file.
  destroyGenEventFile(&file, status);
  if (EXIT_SUCCESS!=*status) return(file);
  file->fptr=NULL;

  // Re-open the file.
  file = openGenEventFile(filename, READWRITE, status);
  if (EXIT_SUCCESS!=*status) return(file);
  
  return(file);
}



GenEventFile* openGenEventFile(const char* const filename,
			       const int mode, int* const status)
{
  GenEventFile* file = newGenEventFile(status);
  if (EXIT_SUCCESS!=*status) return(file);

  headas_chat(4, "open event list file '%s' ...\n", filename);
  if (fits_open_table(&file->fptr, filename, mode, status)) return(file);

  // Determine the row numbers.
  file->row=0;
  fits_get_num_rows(file->fptr, &file->nrows, status);

  // Determine the column numbers.
  if(fits_get_colnum(file->fptr, CASEINSEN, "TIME", &file->ctime, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "PHA", &file->cpha, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "CHARGE", &file->ccharge, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "RAWX", &file->crawx, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "RAWY", &file->crawy, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "FRAME", &file->cframe, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "PILEUP", &file->cpileup, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "PAT_TYPE", &file->cpat_type, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "PAT_ID", &file->cpat_id, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "PAT_ALIG", &file->cpat_alig, status)) 
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
  if (NULL==file->fptr) {
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
  if (fits_write_col(file->fptr, TFLOAT, file->ccharge, file->row, 
		     1, 1, &event->charge, status)) return;
  if (fits_write_col(file->fptr, TINT, file->crawx, file->row, 
		     1, 1, &event->rawx, status)) return;
  if (fits_write_col(file->fptr, TINT, file->crawy, file->row, 
		     1, 1, &event->rawy, status)) return;
  if (fits_write_col(file->fptr, TLONG, file->cframe, file->row, 
		     1, 1, &event->frame, status)) return;
  if (fits_write_col(file->fptr, TINT, file->cpileup, file->row, 
		     1, 1, &event->pileup, status)) return;
  if (fits_write_col(file->fptr, TINT, file->cpat_type, file->row, 
		     1, 1, &event->pat_type, status)) return;
  if (fits_write_col(file->fptr, TINT, file->cpat_id, file->row, 
		     1, 1, &event->pat_id, status)) return;
  if (fits_write_col(file->fptr, TINT, file->cpat_alig, file->row, 
		     1, 1, &event->pat_alig, status)) return;
}



void getGenEventFromFile(const GenEventFile* const file,
			 const int row, GenEvent* const event,
			 int* const status)
{
  // Check if the file has been opened.
  if (NULL==file) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: no GenEventFile opened!\n", *status);
    return;
  }
  if (NULL==file->fptr) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: no GenEventFile opened!\n", *status);
    return;
  }

  // Check if there is still a row available.
  if (row > file->nrows) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: GenEventFile contains no further entries!\n", *status);
    return;
  }

  // Read in the data.
  int anynul = 0;
  event->time = 0.;
  if (0<file->ctime) 
    if (fits_read_col(file->fptr, TDOUBLE, file->ctime, row, 1, 1, 
		      &event->time, &event->time, &anynul, status)) return;
  event->pha = 0;
  if (0<file->cpha) 
    if (fits_read_col(file->fptr, TLONG, file->cpha, row, 1, 1, 
		      &event->pha, &event->pha, &anynul, status)) return;
  event->charge = 0.;
  if (0<file->ccharge) 
    if (fits_read_col(file->fptr, TFLOAT, file->ccharge, row, 1, 1, 
		      &event->charge, &event->charge, &anynul, status)) return;
  event->rawx = 0;
  if (0<file->crawx) 
    if (fits_read_col(file->fptr, TINT, file->crawx, row, 1, 1, 
		      &event->rawx, &event->rawx, &anynul, status)) return;
  event->rawy = 0;
  if (0<file->crawy) 
    if (fits_read_col(file->fptr, TINT, file->crawy, row, 1, 1, 
		      &event->rawy, &event->rawy, &anynul, status)) return;
  event->frame = 0;
  if (0<file->cframe) 
    if (fits_read_col(file->fptr, TLONG, file->cframe, row, 1, 1, 
		      &event->frame, &event->frame, &anynul, status)) return;
  event->pileup = 0;
  if (0<file->cpileup) 
    if (fits_read_col(file->fptr, TINT, file->cpileup, row, 1, 1, 
		      &event->pileup, &event->pileup, &anynul, status)) return;
  event->pat_type = 0;
  if (0<file->cpat_type) 
    if (fits_read_col(file->fptr, TINT, file->cpat_type, row, 1, 1, 
		      &event->pat_type, &event->pat_type, &anynul, status)) return;
  event->pat_id = 0;
  if (0<file->cpat_id) 
    if (fits_read_col(file->fptr, TINT, file->cpat_id, row, 1, 1, 
		      &event->pat_id, &event->pat_id, &anynul, status)) return;
  event->pat_alig = 0;
  if (0<file->cpat_alig) 
    if (fits_read_col(file->fptr, TINT, file->cpat_alig, row, 1, 1, 
		      &event->pat_alig, &event->pat_alig, &anynul, status)) return;

  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: reading from ImpactListFile failed!\n", *status);
    return;
  }

  return;
}



void updateGenEventFromFile(const GenEventFile* const file,
			    const int row, GenEvent* const event,
			    int* const status)
{
  if (fits_write_col(file->fptr, TDOUBLE, file->ctime, row, 
		     1, 1, &event->time, status)) return;
  if (fits_write_col(file->fptr, TLONG, file->cpha, row, 
		     1, 1, &event->pha, status)) return;
  if (fits_write_col(file->fptr, TFLOAT, file->ccharge, row, 
		     1, 1, &event->charge, status)) return;
  if (fits_write_col(file->fptr, TINT, file->crawx, row, 
		     1, 1, &event->rawx, status)) return;
  if (fits_write_col(file->fptr, TINT, file->crawy, row, 
		     1, 1, &event->rawy, status)) return;
  if (fits_write_col(file->fptr, TLONG, file->cframe, row, 
		     1, 1, &event->frame, status)) return;
  if (fits_write_col(file->fptr, TINT, file->cpileup, row, 
		     1, 1, &event->pileup, status)) return;
  if (fits_write_col(file->fptr, TINT, file->cpat_type, row, 
		     1, 1, &event->pat_type, status)) return;
  if (fits_write_col(file->fptr, TINT, file->cpat_id, row, 
		     1, 1, &event->pat_id, status)) return;
  if (fits_write_col(file->fptr, TINT, file->cpat_alig, row, 
		     1, 1, &event->pat_alig, status)) return;
}

