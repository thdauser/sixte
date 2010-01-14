#include "xmseventfile.h"


int openXMSEventFile(XMSEventFile* xef, char* filename, int access_mode)
{
  int status = EXIT_SUCCESS;

  // Call the corresponding routine of the underlying structure.
  status = openEventFile(&xef->generic, filename, access_mode);
  if (EXIT_SUCCESS!=status) return(status);

  // Determine the XMS-specific elements of the event list.
  // Determine the individual column numbers:
  // REQUIRED columns:
  if(fits_get_colnum(xef->generic.fptr, CASEINSEN, "TIME", &xef->ctime, &status)) 
    return(status);
  if(fits_get_colnum(xef->generic.fptr, CASEINSEN, "PHA", &xef->cpha, &status)) 
    return(status);
  if(fits_get_colnum(xef->generic.fptr, CASEINSEN, "RAWX", &xef->crawx, &status)) 
    return(status);
  if(fits_get_colnum(xef->generic.fptr, CASEINSEN, "RAWY", &xef->crawy, &status)) 
    return(status);
  if(fits_get_colnum(xef->generic.fptr, CASEINSEN, "GRADE", &xef->cgrade, &status)) 
    return(status);
  if(fits_get_colnum(xef->generic.fptr, CASEINSEN, "ARRAY", &xef->carray, &status)) 
    return(status);

  return(status);
}



int openNewXMSEventFile(XMSEventFile* xef, char* filename, char* template)
{
  int status=EXIT_SUCCESS;

  // Set the FITS file pointer to NULL. In case that an error occurs during the file
  // generation, we want to avoid that the file pointer points somewhere.
  xef->generic.fptr = NULL;

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
  status = openXMSEventFile(xef, filename, READWRITE);

  return(status);
}



int closeXMSEventFile(XMSEventFile* xef)
{
  // Call the corresponding routine of the underlying structure.
  return(closeEventFile(&xef->generic));
}



int addXMSEvent2File(XMSEventFile* xef, XMSEvent* event)
{
  int status=EXIT_SUCCESS;

  // Insert a new, empty row to the table:
  if (fits_insert_rows(xef->generic.fptr, xef->generic.row, 1, &status)) return(status);
  xef->generic.row++;
  xef->generic.nrows++;

  // Write the event data to the newly created row.
  status=XMSEventFile_writeRow(xef, event, xef->generic.row);

  return(status);
}



int XMSEventFile_getNextRow(XMSEventFile* ef, XMSEvent* event)
{
  int status=EXIT_SUCCESS;

  // Move counter to next line.
  ef->generic.row++;

  // Check if there is still a row available.
  if (ef->generic.row > ef->generic.nrows) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: event list file contains no further entries!\n", status);
    return(status);
  }

  // Read the new XMSEvent from the file.
  status=XMSEventFile_getRow(ef, event, ef->generic.row);

  return(status);
}



int XMSEventFile_getRow(XMSEventFile* ef, XMSEvent* event, long row)
{
  int status=EXIT_SUCCESS;
  int anynul = 0;

  // Check if the specified row is valid.
  if (0==EventFileRowIsValid(&ef->generic, row)) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: invalid row in event file!\n", status);
    return(status);
  }

  // Read in the data.
  event->time = 0.;
  if (fits_read_col(ef->generic.fptr, TDOUBLE, ef->ctime, row, 1, 1, 
		    &event->time, &event->time, &anynul, &status)) return(status);
  event->pha = 0;
  if (fits_read_col(ef->generic.fptr, TLONG, ef->cpha, row, 1, 1, 
		    &event->pha, &event->pha, &anynul, &status)) return(status);
  event->xi = 0;
  if (fits_read_col(ef->generic.fptr, TINT, ef->crawx, row, 1, 1, 
		    &event->xi, &event->xi, &anynul, &status)) return(status);
  event->yi = 0;
  if (fits_read_col(ef->generic.fptr, TINT, ef->crawy, row, 1, 1, 
		    &event->yi, &event->yi, &anynul, &status)) return(status);
  event->grade = 0;
  if (fits_read_col(ef->generic.fptr, TINT, ef->cgrade, row, 1, 1, 
		    &event->grade, &event->grade, &anynul, &status)) return(status);
  event->array = 0;
  if (fits_read_col(ef->generic.fptr, TINT, ef->carray, row, 1, 1, 
		    &event->array, &event->array, &anynul, &status)) return(status);
  
  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: reading from event list failed!\n", status);
    return(status);
  }

  return(status);
}



int XMSEventFile_writeRow(XMSEventFile* xef, XMSEvent* event, long row) {
  int status=EXIT_SUCCESS;

  if (fits_write_col(xef->generic.fptr, TDOUBLE, xef->ctime, row, 
		     1, 1, &event->time, &status)) return(status);
  if (fits_write_col(xef->generic.fptr, TLONG, xef->cpha, row, 
		     1, 1, &event->pha, &status)) return(status);
  if (fits_write_col(xef->generic.fptr, TINT, xef->crawx, row, 
		     1, 1, &event->xi, &status)) return(status);
  if (fits_write_col(xef->generic.fptr, TINT, xef->crawy, row, 
		     1, 1, &event->yi, &status)) return(status);
  if (fits_write_col(xef->generic.fptr, TINT, xef->cgrade, row, 
		     1, 1, &event->grade, &status)) return(status);
  if (fits_write_col(xef->generic.fptr, TINT, xef->carray, row, 
		     1, 1, &event->array, &status)) return(status);

  return(status);
}

