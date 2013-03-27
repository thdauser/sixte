#include "comaeventfile.h"


CoMaEventFile* openCoMaEventFile(char* const filename, 
				 const int access_mode, 
				 int* const status)
{
  //Memory-allocation
  CoMaEventFile* ef=(CoMaEventFile*)malloc(sizeof(CoMaEventFile));
  if (NULL==ef) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("could not allocate memory for CoMaEventFile");
    return(ef);
  }

  //Call the corresponding routine of the underlying structure.
  *status=openEventFile(&ef->generic, filename, access_mode);
  CHECK_STATUS_RET(*status, ef);

  //Determine the CoMa-specific elements of the event list.
  //Determine the individual column numbers:
  //REQUIRED columns:
  fits_get_colnum(ef->generic.fptr, CASEINSEN, "TIME", &ef->ctime, status);
  fits_get_colnum(ef->generic.fptr, CASEINSEN, "CHARGE", &ef->ccharge, status);
  fits_get_colnum(ef->generic.fptr, CASEINSEN, "RAWX", &ef->crawx, status); 
  fits_get_colnum(ef->generic.fptr, CASEINSEN, "RAWY", &ef->crawy, status);
  CHECK_STATUS_RET(*status, ef);

  return(ef);
}


CoMaEventFile* openNewCoMaEventFile(char* const filename, 
				    char* const template, 
				    int* const status)
{
  printf("*** %s \n\n", template);

  // Remove old file if it exists.
  remove(filename);

  // Create a new event list FITS file from a FITS template.
  fitsfile* fptr=NULL;
  char buffer[MAXMSG];
  sprintf(buffer, "%s(%s)", filename, template);
  fits_create_file(&fptr, buffer, status);
  CHECK_STATUS_RET(*status, NULL);

  // Set the time-keyword in the event list header.
  char datestr[MAXMSG];
  int timeref;
  fits_get_system_time(datestr, &timeref, status);
  CHECK_STATUS_RET(*status, NULL);
  fits_update_key(fptr, TSTRING, "DATE", datestr, 
		  "File creation date", status);
  CHECK_STATUS_RET(*status, NULL);

  // Close the newly created file again. It will be immediately re-opened
  // by the standard constructor.
  fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, NULL);

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
  //set internal row counter one element further
  ef->generic.row++;
  //increase number of data-sets in the event-file
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


