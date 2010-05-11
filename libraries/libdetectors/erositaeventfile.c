#include "erositaeventfile.h"


int openeROSITAEventFile(eROSITAEventFile* eef, char* filename, int access_mode)
{
  int status = EXIT_SUCCESS;

  // Call the corresponding routine of the underlying structure.
  status = openEventFile(&eef->generic, filename, access_mode);
  if (EXIT_SUCCESS!=status) return(status);

  // Determine the eROSITA-specific elements of the event list.
  // Determine the individual column numbers:
  // REQUIRED columns:
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "TIME", &eef->ctime, &status)) 
    return(status);
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "PHA", &eef->cpha, &status)) 
    return(status);
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "ENERGY", &eef->cenergy, &status)) 
    return(status);
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "FRAME", &eef->cframe, &status)) 
    return(status);
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "RAWX", &eef->crawx, &status)) 
    return(status);
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "RAWY", &eef->crawy, &status)) 
    return(status);
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "CCDNR", &eef->cccdnr, &status)) 
    return(status);
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "PAT_TYP", &eef->cpat_typ, &status)) 
    return(status);
  if(fits_get_colnum(eef->generic.fptr, CASEINSEN, "PAT_INF", &eef->cpat_inf, &status)) 
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



int openNeweROSITAEventFile(eROSITAEventFile* eef, char* filename, char* template)
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
  } 
  // END of writing time information to Event File FITS header.

  // Add header information about program parameters.
  // The second parameter "1" means that the headers are writte
  // to the first extension.
  HDpar_stamp(fptr, 1, &status);

  // Close the file. It will be re-opened immediately with the
  // standard opening routine.
  if (fits_close_file(fptr, &status)) return(status);


  // Open the newly created FITS file.
  status = openeROSITAEventFile(eef, filename, READWRITE);

  return(status);
}



int closeeROSITAEventFile(eROSITAEventFile* eef)
{
  // Call the corresponding routine of the underlying structure.
  return(closeEventFile(&eef->generic));
}



int addeROSITAEvent2File(eROSITAEventFile* eef, eROSITAEvent* event)
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
  if (fits_write_col(eef->generic.fptr, TFLOAT, eef->cenergy, eef->generic.row, 
		     1, 1, &event->energy, &status)) return(status);
  int rawx = event->xi+1;
  if (fits_write_col(eef->generic.fptr, TINT, eef->crawx, eef->generic.row, 
		     1, 1, &rawx, &status)) return(status);
  int rawy = event->yi+1;
  if (fits_write_col(eef->generic.fptr, TINT, eef->crawy, eef->generic.row, 
		     1, 1, &rawy, &status)) return(status);
  if (fits_write_col(eef->generic.fptr, TLONG, eef->cframe, eef->generic.row, 
		     1, 1, &event->frame, &status)) return(status);
  if (fits_write_col(eef->generic.fptr, TBYTE, eef->cccdnr, eef->generic.row, 
		     1, 1, &event->ccdnr, &status)) return(status);
  if (fits_write_col(eef->generic.fptr, TBYTE, eef->cpat_typ, eef->generic.row, 
		     1, 1, &event->pat_typ, &status)) return(status);
  if (fits_write_col(eef->generic.fptr, TBYTE, eef->cpat_inf, eef->generic.row, 
		     1, 1, &event->pat_inf, &status)) return(status);

  // Set RA, DEC, X, and Y to default values.
  event->ra = SIXT_NAN;
  if (fits_write_col(eef->generic.fptr, TDOUBLE, eef->cra, eef->generic.row, 
		     1, 1, &event->ra, &status)) return(status);
  event->dec = SIXT_NAN;
  if (fits_write_col(eef->generic.fptr, TDOUBLE, eef->cdec, eef->generic.row, 
		     1, 1, &event->dec, &status)) return(status);
  event->sky_xi = 0;
  if (fits_write_col(eef->generic.fptr, TLONG, eef->cskyx, eef->generic.row, 
		     1, 1, &event->sky_xi, &status)) return(status);
  event->sky_yi = 0;
  if (fits_write_col(eef->generic.fptr, TLONG, eef->cskyy, eef->generic.row, 
		     1, 1, &event->sky_yi, &status)) return(status);

  return(status);
}



int eROSITAEventFile_getNextRow(eROSITAEventFile* eef, eROSITAEvent* event)
{
  int status=EXIT_SUCCESS;

  // Move counter to next line.
  eef->generic.row++;

  // Check if there is still a row available.
  if (eef->generic.row > eef->generic.nrows) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: event list file contains no further entries!\n", status);
    return(status);
  }

  // Read the new HTRSEvent from the file.
  status=eROSITAEventFile_getRow(eef, event, eef->generic.row);

  return(status);
}



int eROSITAEventFile_getRow(eROSITAEventFile* eef, eROSITAEvent* event, long row)
{
  int status=EXIT_SUCCESS;
  int anynul = 0;

  // Check if there is still a row available.
  if (row > eef->generic.nrows) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: event list file does not contain the requested line!\n", 
		   status);
    return(status);
  }

  // Read in the data.
  event->time = 0.;
  if (fits_read_col(eef->generic.fptr, TDOUBLE, eef->ctime, row, 1, 1, 
		    &event->time, &event->time, &anynul, &status)) return(status);
  event->pha = 0;
  if (fits_read_col(eef->generic.fptr, TLONG, eef->cpha, row, 1, 1, 
		    &event->pha, &event->pha, &anynul, &status)) return(status);
  event->energy = 0.;
  if (fits_read_col(eef->generic.fptr, TFLOAT, eef->cenergy, row, 1, 1, 
		    &event->energy, &event->energy, &anynul, &status)) return(status);

  event->xi = 0;
  if (fits_read_col(eef->generic.fptr, TINT, eef->crawx, row, 1, 1, 
		    &event->xi, &event->xi, &anynul, &status)) return(status);
  event->xi--;
  event->yi = 0;
  if (fits_read_col(eef->generic.fptr, TINT, eef->crawy, row, 1, 1, 
		    &event->yi, &event->yi, &anynul, &status)) return(status);
  event->yi--;

  event->frame = 0;
  if (fits_read_col(eef->generic.fptr, TLONG, eef->cframe, row, 1, 1, 
		    &event->frame, &event->frame, &anynul, &status)) return(status);
  event->ccdnr = 0;
  if (fits_read_col(eef->generic.fptr, TBYTE, eef->cccdnr, row, 1, 1, 
		    &event->ccdnr, &event->ccdnr, &anynul, &status)) return(status);
  event->pat_typ = 0;
  if (fits_read_col(eef->generic.fptr, TBYTE, eef->cpat_typ, row, 1, 1, 
		    &event->pat_typ, &event->pat_typ, &anynul, &status)) return(status);
  event->pat_inf = 0;
  if (fits_read_col(eef->generic.fptr, TBYTE, eef->cpat_inf, row, 1, 1, 
		    &event->pat_inf, &event->pat_inf, &anynul, &status)) return(status);

  event->ra = SIXT_NAN;
  if (fits_read_col(eef->generic.fptr, TDOUBLE, eef->cra, row, 1, 1, 
		    &event->ra, &event->ra, &anynul, &status)) return(status);
  event->dec = SIXT_NAN;
  if (fits_read_col(eef->generic.fptr, TDOUBLE, eef->cdec, row, 1, 1, 
		    &event->dec, &event->dec, &anynul, &status)) return(status);
  event->sky_xi = 0;
  if (fits_read_col(eef->generic.fptr, TLONG, eef->cskyx, row, 1, 1, 
		    &event->sky_xi, &event->sky_xi, &anynul, &status)) return(status);
  event->sky_yi = 0;
  if (fits_read_col(eef->generic.fptr, TLONG, eef->cskyy, row, 1, 1, 
		    &event->sky_yi, &event->sky_yi, &anynul, &status)) return(status);

  
  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: reading from event list failed!\n", status);
    return(status);
  }

  return(status);
}


