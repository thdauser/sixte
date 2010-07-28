#include "pnccdeventfile.h"


int openpnCCDEventFile(pnCCDEventfile* pef, char* filename, int access_mode)
{
	int status = EXIT_SUCCESS;

	// Open the EventFile
	// TODO: openEventFile has to be written (or maybe not...)
	status = openEventFile(&pef->generic, filename, access_mode);
	if (EXIT_SUCCESS!=status) return(status);

	// Determine the pnCCD-specific elements of the event list
	// Determine the individual column numbers:
	if(fits_get_colnum(pef->generic.fptr, CASEINSEN, "TIME", &pef->ctime, &status)) return(status);
	// !!! NOTE: get the colnums of all the elements of the event list

	//TODO: Determine what parameters are needed for the pnCCD Eventfile

	return(status);
}
	


int openNewpnCCDEventFile(pnCCDEventFile* pef, char* filename, char* template)
{
	int status=EXIT_SUCCESS;

	// Set the FITS file pointer to NULL. In case that an error occurs during the file
  // generation, we want to avoid that the file pointer points somewhere.
  pef->generic.fptr = NULL;

	// Remove old file if it exists
	remove(filename);

	// Create a new event list FITS file from a FITS template
	fitsfile* fptr=NULL;
	char buffer[MAXMSG];
	sprintf(buffer, "%s(%s)", filename, template);
	if (fits_create_file(&fptr, buffer, &status)) return(status);

	// Set the time-keyword in the Event List Header
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
  status = openpnCCDEventFile(pef, filename, READWRITE);

  return(status);
}



int closepnCCDEventFile(pnCCDEventFile* pef)
{
	return(closeEventFile(&pef->generic));
}



//TODO: Check if the pnCCDEvent structure exists
int addpnCCDEvent2File(pnCCDEventFile* pef, pnCCDEvent* event)
{
	int status=EXIT_SUCCESS;

	// Insert a new, empty row to the table:
	if (fits_insert_rows(pef->generic.fptr, pef->generic.row, 1, &status)) return(status);
	pef->generic.row++;
	pef->generic.nrows++;

	if (fits_write_col(pef->generic.fptr, TDOUBLE, pef->ctime, pef->generic.row, 1, 1, &event->time, &status)) return(status);

	//TODO: analog to the get_colnum function above the parameters has to be
	//			defined

	return(status);
}

int pnCCDEventFile_getNextRow(pnCCDEventFile* pef, pnCCDEvent* event)
{
	int status=EXIT_SUCCESS;

	// Move counter to next line
	pef->generic.row++;

	// Check if there is still a row available.
	if (pef->generic.row > pef->generic.nrows) {
		status = EXIT_FAILURE;
		HD_ERROR_THROW("Error: event list file contains no further enries!\n", status);
		return(status);
	}

	// Read the new pnCCDEvent from the file
	status=pnCCDEventFile_getRow(pef, event, pef->generic.row);

	return(status);
}



int pnCCDEventFile_getRow(pnCCDEventFile* pef, pnCCDEvent* event, long row)
{
	int status=EXIT_SUCCESS;
	int anynul = 0;

	// Check if there is still a row available
	if (rwo > pef->generic.nrows) {
		status = EXIT_FAILURE;
		HD_ERROR_THROW("Error: Event list file does not contain the requested line!\n", status);
		return(status);
	}

	// Read in the data 
	event->time = 0.;
	if (fits_read_col(pef->generic.fptr, TDOUBLE, pef->ctime, row, 1, 1,
				&event->time, &event->time, &anynul, &status)) return(status);

	// TODO: Analog to the TODOs above, define parameters

	if (0!=anynul) {
		status = EXIT_FAILURE;
		HD_ERROR_THROW("Error: reading from event list failed!\n", status);
		return(status);
	}

	return(status);
}
