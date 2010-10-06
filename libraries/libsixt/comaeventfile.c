#include "comaeventfile.h"


CoMaEventFile* openCoMaEventFile(char* filename, int access_mode, int* status)
{
  CoMaEventFile* ef=(CoMaEventFile*)malloc(sizeof(CoMaEventFile));
  if (NULL==ef) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for CoMaEventFile "
		   "object!\n", *status);
    return(ef);
  }
  
  // Call the corresponding routine of the underlying structure.
  *status = openEventFile(&ef->generic, filename, access_mode);
  if (EXIT_SUCCESS!=*status) return(ef);

  // Determine the CoMa-specific elements of the event list.
  // Determine the individual column numbers:
  // REQUIRED columns:
  if(fits_get_colnum(ef->generic.fptr, CASEINSEN, "TIME", &ef->ctime, status)) 
    return(ef);
  if(fits_get_colnum(ef->generic.fptr, CASEINSEN, "CHARGE", &ef->ccharge, 
		     status)) 
    return(ef);
  if(fits_get_colnum(ef->generic.fptr, CASEINSEN, "RAWX", &ef->crawx, status)) 
    return(ef);
  if(fits_get_colnum(ef->generic.fptr, CASEINSEN, "RAWY", &ef->crawy, status)) 
    return(ef);

  return(ef);
}



CoMaEventFile* openNewCoMaEventFile(char* filename, char* template, int* status)
{
  // Remove old file if it exists.
  remove(filename);
  // Create a new event list FITS file from a FITS template.
  fitsfile* fptr=NULL;
  char buffer[MAXMSG];
  sprintf(buffer, "%s(%s)", filename, template);
  if (fits_create_file(&fptr, buffer, status)) return(NULL);

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
	  return(NULL);
      }
    }
  } // END of writing time information to Event File FITS header.

  // Close the newly created file again. It will be immediately re-opened
  // by the standard constructor.
  if (fits_close_file(fptr, status)) return(NULL);


  // Open the newly created FITS file.
  return(openCoMaEventFile(filename, READWRITE, status));
}



int closeCoMaEventFile(CoMaEventFile* ef)
{
  // Call the corresponding routine of the underlying structure.
  return(closeEventFile(&ef->generic));
}



int addCoMaEvent2File(CoMaEventFile* ef, CoMaEvent* event)
{
  int status=EXIT_SUCCESS;

  // Insert a new, empty row to the table:
  if (fits_insert_rows(ef->generic.fptr, ef->generic.row, 1, &status)) 
    return(status);
  ef->generic.row++;
  ef->generic.nrows++;

  if (fits_write_col(ef->generic.fptr, TDOUBLE, ef->ctime, ef->generic.row, 
		     1, 1, &event->time, &status)) return(status);
  if (fits_write_col(ef->generic.fptr, TFLOAT, ef->ccharge, ef->generic.row, 
		     1, 1, &event->charge, &status)) return(status);
  if (fits_write_col(ef->generic.fptr, TINT, ef->crawx, ef->generic.row, 
		     1, 1, &event->rawx, &status)) return(status);
  if (fits_write_col(ef->generic.fptr, TINT, ef->crawy, ef->generic.row, 
		     1, 1, &event->rawy, &status)) return(status);

  return(status);
}



int CoMaEventFile_getNextRow(CoMaEventFile* ef, CoMaEvent* event)
{
  int status=EXIT_SUCCESS;
  int anynul=0;

  // Move counter to next line.
  ef->generic.row++;

  // Check if there is still a row available.
  if (ef->generic.row > ef->generic.nrows) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: event list file contains no further entries!\n", 
		   status);
    return(status);
  }

  // Read in the data.
  event->time = 0.;
  if (fits_read_col(ef->generic.fptr, TDOUBLE, ef->ctime, ef->generic.row, 1, 
		    1, &event->time, &event->time, &anynul, &status)) 
    return(status);
  event->charge = 0.;
  if (fits_read_col(ef->generic.fptr, TLONG, ef->ccharge, ef->generic.row, 1, 
		    1, &event->charge, &event->charge, &anynul, &status)) 
    return(status);
  event->rawx = 0;
  if (fits_read_col(ef->generic.fptr, TINT, ef->crawx, ef->generic.row, 1, 1, 
		    &event->rawx, &event->rawx, &anynul, &status)) return(status);
  event->rawy = 0;
  if (fits_read_col(ef->generic.fptr, TINT, ef->crawy, ef->generic.row, 1, 1, 
		    &event->rawy, &event->rawy, &anynul, &status)) return(status);
  
  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: reading from event list failed!\n", status);
    return(status);
  }

  return(status);
}


