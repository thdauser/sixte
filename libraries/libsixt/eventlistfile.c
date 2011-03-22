#include "eventlistfile.h"


EventListFile* newEventListFile(int* const status)
{
  EventListFile* file = (EventListFile*)malloc(sizeof(EventListFile));
  if (NULL==file) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for EventListFile failed!\n", 
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
  file->cph_id =0;
  file->csrc_id=0;

  return(file);
}



void freeEventListFile(EventListFile** const file, int* const status)
{
  if (NULL!=*file) {
    if (NULL!=(*file)->fptr) {
      fits_close_file((*file)->fptr, status);
    }
    free(*file);
    *file=NULL;
  }
}



EventListFile* openNewEventListFile(const char* const filename,
				    const char* const template,
				    int* const status)
{
  EventListFile* file = newEventListFile(status);
  if (EXIT_SUCCESS!=*status) return(file);

  // Remove old file, if it exists.
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
  freeEventListFile(&file, status);
  if (EXIT_SUCCESS!=*status) return(file);

  // Re-open the file.
  file = openEventListFile(filename, READWRITE, status);
  if (EXIT_SUCCESS!=*status) return(file);
  
  return(file);
}



EventListFile* openEventListFile(const char* const filename,
				 const int mode, int* const status)
{
  EventListFile* file = newEventListFile(status);
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
  if(fits_get_colnum(file->fptr, CASEINSEN, "PH_ID", &file->cph_id, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "SRC_ID", &file->csrc_id, status)) 
    return(file);

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
    *status = EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "Error: the maximum number of photons contributing "
	    "to a single event is different for the parameter set "
	    "in the simulation (%d) and in the event list "
	    "template file (%ld)!\n", NEVENTPHOTONS, repeat);
    HD_ERROR_THROW(msg, *status);
    return(file);
  }
  // SRC_ID.
  fits_get_coltype(file->fptr, file->csrc_id, &typecode, &repeat,
		   &width, status);
  CHECK_STATUS_RET(*status, file);
  if (repeat!=NEVENTPHOTONS) {
    // Throw an error.
    *status = EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "Error: the maximum number of photons contributing "
	    "to a single event is different for the parameter set "
	    "in the simulation (%d) and in the event list "
	    "template file (%ld)!\n", NEVENTPHOTONS, repeat);
    HD_ERROR_THROW(msg, *status);
    return(file);
  }

  return(file);
}



void addEvent2File(EventListFile* const file, Event* const event, 
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
  if (fits_write_col(file->fptr, TLONG, file->cph_id, file->row, 
		     1, NEVENTPHOTONS, &event->ph_id, status)) return;
  if (fits_write_col(file->fptr, TLONG, file->csrc_id, file->row, 
		     1, NEVENTPHOTONS, &event->src_id, status)) return;
}



void getEventFromFile(const EventListFile* const file,
		      const int row, Event* const event,
		      int* const status)
{
  // Check if the file has been opened.
  if (NULL==file) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: no EventListFile opened!\n", *status);
    return;
  }
  if (NULL==file->fptr) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: no EventListFile opened!\n", *status);
    return;
  }

  // Check if there is still a row available.
  if (row > file->nrows) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: EventListFile contains no further entries!\n", *status);
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
  CHECK_STATUS_VOID(*status);

  fits_read_col(file->fptr, TLONG, file->cpha, row, 1, 1, 
		&lnull, &event->pha, &anynul, status);
  CHECK_STATUS_VOID(*status);

  fits_read_col(file->fptr, TFLOAT, file->ccharge, row, 1, 1, 
		&fnull, &event->charge, &anynul, status);
  CHECK_STATUS_VOID(*status);

  fits_read_col(file->fptr, TINT, file->crawx, row, 1, 1, 
		&inull, &event->rawx, &anynul, status);
  CHECK_STATUS_VOID(*status);

  fits_read_col(file->fptr, TINT, file->crawy, row, 1, 1, 
		&inull, &event->rawy, &anynul, status);
  CHECK_STATUS_VOID(*status);

  fits_read_col(file->fptr, TLONG, file->cframe, row, 1, 1, 
		&lnull, &event->frame, &anynul, status);
  CHECK_STATUS_VOID(*status);

  fits_read_col(file->fptr, TLONG, file->cph_id, row, 1, NEVENTPHOTONS, 
		&lnull, &event->ph_id, &anynul, status);
  CHECK_STATUS_VOID(*status);

  fits_read_col(file->fptr, TLONG, file->csrc_id, row, 1, NEVENTPHOTONS, 
		&lnull, &event->src_id, &anynul, status);
  CHECK_STATUS_VOID(*status);

  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: reading from ImpactListFile failed!\n", *status);
    return;
  }

  return;
}



void updateEventInFile(const EventListFile* const file,
		       const int row, Event* const event,
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
}

