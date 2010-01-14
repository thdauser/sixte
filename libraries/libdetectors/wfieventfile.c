#include "wfieventfile.h"


int openWFIEventFile(WFIEventFile* wef, char* filename, int access_mode)
{
  int status = EXIT_SUCCESS;

  // Call the corresponding routine of the underlying structure.
  status = openEventFile(&wef->generic, filename, access_mode);
  if (EXIT_SUCCESS!=status) return(status);

  // Determine the WFI-specific elements of the event list.
  // Determine the individual column numbers:
  // REQUIRED columns:
  if(fits_get_colnum(wef->generic.fptr, CASEINSEN, "TIME", &wef->ctime, &status)) 
    return(status);
  if(fits_get_colnum(wef->generic.fptr, CASEINSEN, "PHA", &wef->cpha, &status)) 
    return(status);
  if(fits_get_colnum(wef->generic.fptr, CASEINSEN, "FRAME", &wef->cframe, &status)) 
    return(status);
  if(fits_get_colnum(wef->generic.fptr, CASEINSEN, "COLUMN", &wef->crawx, &status)) 
    return(status);
  if(fits_get_colnum(wef->generic.fptr, CASEINSEN, "ROW", &wef->crawy, &status)) 
    return(status);

  if(fits_get_colnum(wef->generic.fptr, CASEINSEN, "PATNUM", &wef->cpatnum, &status)) 
    return(status);
  if(fits_get_colnum(wef->generic.fptr, CASEINSEN, "PATID", &wef->cpatid, &status)) 
    return(status);
  if(fits_get_colnum(wef->generic.fptr, CASEINSEN, "PILEUP", &wef->cpileup, &status)) 
    return(status);

  // Read the header keywords:
  char comment[MAXMSG];
  if(fits_read_key(wef->generic.fptr, TLONG, "COLUMNS", &wef->columns, comment, &status))
    return(status);
  if(fits_read_key(wef->generic.fptr, TLONG, "ROWS", &wef->rows, comment, &status))
    return(status);

  return(status);
}



int openNewWFIEventFile(WFIEventFile* wef, char* filename, char* template)
{
  int status=EXIT_SUCCESS;

  // Set the FITS file pointer to NULL. In case that an error occurs during the file
  // generation, we want to avoid that the file pointer points somewhere.
  wef->generic.fptr = NULL;

  // Remove old file if it exists.
  remove(filename);

  // Create a new event list FITS file from a FITS template.
  fitsfile* fptr=NULL;
  char buffer[MAXMSG];
  sprintf(buffer, "%s(%s)", filename, template);
  if (fits_create_file(&fptr, buffer, &status)) return(status);

  // Set the time-keyword in the Event List Header.
  // See also: Stevens, "Advanced Programming in the UNIX environment", p. 155 ff.
  time_t current_time;
  if (0 != time(&current_time)) {
    struct tm* current_time_utc = gmtime(&current_time);
    if (NULL != current_time_utc) {
      char current_time_str[MAXMSG];
      if (strftime(current_time_str, MAXMSG, "%Y-%m-%dT%H:%M:%S", current_time_utc) > 0) {
	// Return value should be == 19 !
	if (fits_update_key(fptr, TSTRING, "DATE-OBS", current_time_str, 
			    "Start Time (UTC) of exposure", &status)) return(status);
      }
    }
  } // END of writing time information to Event File FITS header.

  // Add header information about program parameters.
  // The second parameter "1" means that the headers are writte
  // to the first extension.
  HDpar_stamp(fptr, 1, &status);

  // Close the file. It will be re-opened immediately with the
  // standard opening routine.
  if (fits_close_file(fptr, &status)) return(status);


  // Open the newly created FITS file.
  status = openWFIEventFile(wef, filename, READWRITE);

  return(status);
}



int closeWFIEventFile(WFIEventFile* wef)
{
  // Call the corresponding routine of the underlying structure.
  return(closeEventFile(&wef->generic));
}



int addWFIEvent2File(WFIEventFile* wef, WFIEvent* event)
{
  int status=EXIT_SUCCESS;

  // Insert a new, empty row to the table:
  if (fits_insert_rows(wef->generic.fptr, wef->generic.row, 1, &status)) return(status);
  wef->generic.row++;
  wef->generic.nrows++;

  if (fits_write_col(wef->generic.fptr, TDOUBLE, wef->ctime, wef->generic.row, 
		     1, 1, &event->time, &status)) return(status);
  if (fits_write_col(wef->generic.fptr, TLONG, wef->cpha, wef->generic.row, 
		     1, 1, &event->pha, &status)) return(status);
  if (fits_write_col(wef->generic.fptr, TINT, wef->crawx, wef->generic.row, 
		     1, 1, &event->xi, &status)) return(status);
  if (fits_write_col(wef->generic.fptr, TINT, wef->crawy, wef->generic.row, 
		     1, 1, &event->yi, &status)) return(status);
  if (fits_write_col(wef->generic.fptr, TLONG, wef->cframe, wef->generic.row, 
		     1, 1, &event->frame, &status)) return(status);

  if (fits_write_col(wef->generic.fptr, TLONG, wef->cpatnum, wef->generic.row, 
		     1, 1, &event->patnum, &status)) return(status);
  if (fits_write_col(wef->generic.fptr, TLONG, wef->cpatid, wef->generic.row, 
		     1, 1, &event->patid, &status)) return(status);
  if (fits_write_col(wef->generic.fptr, TLONG, wef->cpileup, wef->generic.row, 
		     1, 1, &event->pileup, &status)) return(status);

  return(status);
}



int WFIEventFile_getNextRow(WFIEventFile* ef, WFIEvent* event)
{
  int status=EXIT_SUCCESS;
  int anynul = 0;

  // Move counter to next line.
  ef->generic.row++;

  // Check if there is still a row available.
  if (ef->generic.row > ef->generic.nrows) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: event list file contains no further entries!\n", status);
    return(status);
  }

  // Read in the data.
  event->time = 0.;
  if (fits_read_col(ef->generic.fptr, TDOUBLE, ef->ctime, ef->generic.row, 1, 1, 
		    &event->time, &event->time, &anynul, &status)) return(status);
  event->pha = 0;
  if (fits_read_col(ef->generic.fptr, TLONG, ef->cpha, ef->generic.row, 1, 1, 
		    &event->pha, &event->pha, &anynul, &status)) return(status);
  event->xi = 0;
  if (fits_read_col(ef->generic.fptr, TINT, ef->crawx, ef->generic.row, 1, 1, 
		    &event->xi, &event->xi, &anynul, &status)) return(status);
  event->yi = 0;
  if (fits_read_col(ef->generic.fptr, TINT, ef->crawy, ef->generic.row, 1, 1, 
		    &event->yi, &event->yi, &anynul, &status)) return(status);
  event->frame = 0;
  if (fits_read_col(ef->generic.fptr, TLONG, ef->cframe, ef->generic.row, 1, 1, 
		    &event->frame, &event->frame, &anynul, &status)) return(status);
  event->patnum = 0;
  if (fits_read_col(ef->generic.fptr, TLONG, ef->cpatnum, ef->generic.row, 1, 1, 
		    &event->patnum, &event->patnum, &anynul, &status)) return(status);
  event->patid = 0;
  if (fits_read_col(ef->generic.fptr, TLONG, ef->cpatid, ef->generic.row, 1, 1, 
		    &event->patid, &event->patid, &anynul, &status)) return(status);
  event->pileup = 0;
  if (fits_read_col(ef->generic.fptr, TLONG, ef->cpileup, ef->generic.row, 1, 1, 
		    &event->pileup, &event->pileup, &anynul, &status)) return(status);
  
  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: reading from event list failed!\n", status);
    return(status);
  }

  return(status);
}


