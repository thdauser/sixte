#include "htrseventfile.h"


int openHTRSEventFile(HTRSEventFile* hef, char* filename, int access_mode)
{
  int status = EXIT_SUCCESS;

  // Call the corresponding routine of the underlying structure.
  status = openEventFile(&hef->generic, filename, access_mode);
  if (EXIT_SUCCESS!=status) return(status);

  // Determine the HTRS-specific elements of the event list.
  // Determine the individual column numbers:
  // REQUIRED columns:
  if(fits_get_colnum(hef->generic.fptr, CASEINSEN, "TIME", &hef->ctime, &status)) 
    return(status);
  if(fits_get_colnum(hef->generic.fptr, CASEINSEN, "PHA", &hef->cpha, &status)) 
    return(status);
  if(fits_get_colnum(hef->generic.fptr, CASEINSEN, "ENERGY", &hef->cenergy, &status)) 
    return(status);
  if(fits_get_colnum(hef->generic.fptr, CASEINSEN, "PIXEL", &hef->cpixel, &status)) 
    return(status);
  if(fits_get_colnum(hef->generic.fptr, CASEINSEN, "GRADE", &hef->cgrade, &status)) 
    return(status);
  if(fits_get_colnum(hef->generic.fptr, CASEINSEN, "X", &hef->cx, &status)) 
    return(status);
  if(fits_get_colnum(hef->generic.fptr, CASEINSEN, "Y", &hef->cy, &status)) 
    return(status);

  return(status);
}



int openNewHTRSEventFile(HTRSEventFile* hef, char* filename, char* template)
{
  int status=EXIT_SUCCESS;

  // Set the FITS file pointer to NULL. In case that an error occurs during the file
  // generation, we want to avoid that the file pointer points somewhere.
  hef->generic.fptr = NULL;

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
  status = openHTRSEventFile(hef, filename, READWRITE);

  return(status);
}



int closeHTRSEventFile(HTRSEventFile* hef)
{
  // Call the corresponding routine of the underlying structure.
  return(closeEventFile(&hef->generic));
}



int addHTRSEvent2File(HTRSEventFile* hef, HTRSEvent* event)
{
  int status=EXIT_SUCCESS;

  // Insert a new, empty row to the table:
  if (fits_insert_rows(hef->generic.fptr, hef->generic.row, 1, &status)) return(status);
  hef->generic.row++;
  hef->generic.nrows++;

  if (fits_write_col(hef->generic.fptr, TDOUBLE, hef->ctime, hef->generic.row, 
		     1, 1, &event->time, &status)) return(status);
  if (fits_write_col(hef->generic.fptr, TLONG, hef->cpha, hef->generic.row, 
		     1, 1, &event->pha, &status)) return(status);
  if (fits_write_col(hef->generic.fptr, TFLOAT, hef->cenergy, hef->generic.row, 
		     1, 1, &event->energy, &status)) return(status);
  int pixel = event->pixel+1;
  if (fits_write_col(hef->generic.fptr, TINT, hef->cpixel, hef->generic.row, 
		     1, 1, &pixel, &status)) return(status);

  if (fits_write_col(hef->generic.fptr, TINT, hef->cgrade, hef->generic.row, 
		     1, 1, &event->grade, &status)) return(status);
  if (fits_write_col(hef->generic.fptr, TDOUBLE, hef->cx, hef->generic.row, 
		     1, 1, &event->x, &status)) return(status);
  if (fits_write_col(hef->generic.fptr, TDOUBLE, hef->cy, hef->generic.row, 
		     1, 1, &event->y, &status)) return(status);
  return(status);
}



int HTRSEventFile_getNextRow(HTRSEventFile* hef, HTRSEvent* event)
{
  int status=EXIT_SUCCESS;

  // Move counter to next line.
  hef->generic.row++;

  // Check if there is still a row available.
  if (hef->generic.row > hef->generic.nrows) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: event list file contains no further entries!\n", status);
    return(status);
  }

  // Read the new HTRSEvent from the file.
  status=HTRSEventFile_getRow(hef, event, hef->generic.row);

  return(status);
}



int HTRSEventFile_getRow(HTRSEventFile* hef, HTRSEvent* event, long row)
{
  int status=EXIT_SUCCESS;
  int anynul=0;

  // Check if there is still a row available.
  if (row > hef->generic.nrows) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: event list file does not contain the requested line!\n", status);
    return(status);
  }

  // Read in the data.
  event->time = 0.;
  if (fits_read_col(hef->generic.fptr, TDOUBLE, hef->ctime, row, 1, 1, 
		    &event->time, &event->time, &anynul, &status)) return(status);
  event->pha = 0;
  if (fits_read_col(hef->generic.fptr, TLONG, hef->cpha, row, 1, 1, 
		    &event->pha, &event->pha, &anynul, &status)) return(status);
  event->energy = 0;
  if (fits_read_col(hef->generic.fptr, TFLOAT, hef->cenergy, row, 1, 1, 
		    &event->energy, &event->energy, &anynul, &status)) return(status);
  event->pixel = 0;
  if (fits_read_col(hef->generic.fptr, TINT, hef->cpixel, row, 1, 1, 
		    &event->pixel, &event->pixel, &anynul, &status)) return(status);
  event->pixel--;

  event->grade = 0;
  if (fits_read_col(hef->generic.fptr, TINT, hef->cgrade, row, 1, 1, 
		    &event->grade, &event->grade, &anynul, &status)) return(status);
  event->x = 0.;
  if (fits_read_col(hef->generic.fptr, TDOUBLE, hef->cx, row, 1, 1, 
		    &event->x, &event->x, &anynul, &status)) return(status);
  event->y = 0.;
  if (fits_read_col(hef->generic.fptr, TDOUBLE, hef->cy, row, 1, 1, 
		    &event->y, &event->y, &anynul, &status)) return(status);
  
  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: reading from event list failed!\n", status);
    return(status);
  }

  return(status);
}

