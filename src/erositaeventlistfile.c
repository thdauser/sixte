#include "erositaeventlistfile.h"



int openeROSITAEventlistFile(eROSITAEventlistFile* eef, char* filename, int access_mode)
{
  int status = EXIT_SUCCESS;

  // Call the corresponding routine of the underlying structure.
  status = openEventlistFile(&eef->generic, filename, access_mode);
  if (EXIT_SUCCESS!=status) return(status);

  // Determine the eROSITA-specific elements of the event list.
  // Determine the individual column numbers:
  // REQUIRED columns:
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "TIME", &eef->ctime, &status)) 
    return(status);
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "PHA", &eef->cpha, &status)) 
    return(status);
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "FRAME", &eef->cframe, &status)) 
    return(status);
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "RAWX", &eef->crawx, &status)) 
    return(status);
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "RAWY", &eef->crawy, &status)) 
    return(status);

  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "RA", &eef->cra, &status)) 
    return(status);
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "DEC", &eef->cdec, &status)) 
    return(status);
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "X", &eef->cskyx, &status)) 
    return(status);
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "Y", &eef->cskyy, &status)) 
    return(status);

  return(status);
}


int openNeweROSITAEventlistFile(eROSITAEventlistFile* eef, char* filename, char* template)
{
  int status=EXIT_SUCCESS;

  // Set the FITS file pointer to NULL. In case that an error occurs during the file
  // generation, we want to avoid that the file pointer points somewhere.
  eef->generic.fptr = NULL;

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
  status = openeROSITAEventlistFile(eef, filename, READWRITE);

  return(status);
}



int closeeROSITAEventlistFile(eROSITAEventlistFile* eef)
{
  // Call the corresponding routine of the underlying structure.
  return(closeEventlistFile(&eef->generic));
}



int addeROSITAEvent2File(eROSITAEventlistFile* eef, eROSITAEvent* event)
{
  int status=EXIT_SUCCESS;

  // Insert a new, empty row to the table:
  if (fits_insert_rows(eef->generic.fptr, eef->generic.row, 1, &status)) return(status);
  eef->generic.row++;
  eef->generic.nrows++;

  if (fits_write_col(eef->generic.fptr, TDOUBLE, eef->ctime, eef->generic.row, 
		     1, 1, &event->time, &status)) return(status);
  if (fits_write_col(eef->generic.fptr, TLONG, eef->cpha, eef->generic.row, 
		     1, 1, &event->pha, &status)) return(status);
  if (fits_write_col(eef->generic.fptr, TINT, eef->crawx, eef->generic.row, 
		     1, 1, &event->xi, &status)) return(status);
  if (fits_write_col(eef->generic.fptr, TINT, eef->crawy, eef->generic.row, 
		     1, 1, &event->yi, &status)) return(status);
  if (fits_write_col(eef->generic.fptr, TLONG, eef->cframe, eef->generic.row, 
		     1, 1, &event->frame, &status)) return(status);

  // Set RA, DEC, X, and Y to default values.
  event->ra = NAN;
  if (fits_write_col(eef->generic.fptr, TLONG, eef->cra, eef->generic.row, 
		     1, 1, &event->ra, &status)) return(status);
  event->dec = NAN;
  if (fits_write_col(eef->generic.fptr, TLONG, eef->cdec, eef->generic.row, 
		     1, 1, &event->dec, &status)) return(status);
  event->sky_xi = 0;
  if (fits_write_col(eef->generic.fptr, TLONG, eef->cskyx, eef->generic.row, 
		     1, 1, &event->sky_xi, &status)) return(status);
  event->sky_yi = 0;
  if (fits_write_col(eef->generic.fptr, TLONG, eef->cskyy, eef->generic.row, 
		     1, 1, &event->sky_yi, &status)) return(status);

  return(status);
}


