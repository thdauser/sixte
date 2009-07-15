#include "wfieventlistfile.h"



int openWFIEventlistFile(WFIEventlistFile* wef, char* filename, int access_mode)
{
  int status = EXIT_SUCCESS;

  // Call the corresponding routine of the underlying structure.
  status = openEventlistFile(&wef->generic, filename, access_mode);
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

  return(status);
}


int openNewWFIEventlistFile(WFIEventlistFile* wef, char* filename, char* template)
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

  if (fits_close_file(fptr, &status)) return(status);


  // Open the newly created FITS file.
  status = openWFIEventlistFile(wef, filename, READWRITE);

  return(status);
}



int closeWFIEventlistFile(WFIEventlistFile* wef)
{
  // Call the corresponding routine of the underlying structure.
  return(closeEventlistFile(&wef->generic));
}


int addWFIEvent2File(WFIEventlistFile* wef, WFIEvent* event)
{
  int status=EXIT_SUCCESS;

  // Insert a new, empty row to the table:
  if (fits_insert_rows(wef->generic.fptr, wef->generic.row, 1, &status)) {
    if (EXIT_SUCCESS!=status) {
      printf("remove this!\n");
      return(status);
    } else {
      printf("change this!\n");
    }
  }
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

  // Set default values for PATNUM and PATID:
  // PATID has to be set to -1 !! 
  // Otherwise the pattern recognition algorithm doesn't work properly.
  // PATNUM
  event->patnum = 0;
  if (fits_write_col(wef->generic.fptr, TLONG, wef->cpatnum, wef->generic.row, 
		     1, 1, &event->patnum, &status)) return(status);
  // PATID
  event->patid = -1;
  if (fits_write_col(wef->generic.fptr, TLONG, wef->cpatid, wef->generic.row, 
		     1, 1, &event->patid, &status)) return(status);
  // Pile-up
  event->pileup = 0;
  if (fits_write_col(wef->generic.fptr, TLONG, wef->cpileup, wef->generic.row, 
		     1, 1, &event->pileup, &status)) return(status);

  return(status);
}


