#include "ladeventlistfile.h"


LADEventListFile* newLADEventListFile(int* const status)
{
  LADEventListFile* file=
    (LADEventListFile*)malloc(sizeof(LADEventListFile));
  CHECK_NULL_RET(file, *status, 
		 "memory allocation for LADEventListFile failed", file);

  // Initialize pointers with NULL.
  file->fptr=NULL;

  // Initialize values.
  file->nrows=0;
  file->row  =0;
  file->ctime=0;
  file->csignal=0;
  file->cpanel =0;
  file->cmodule=0;
  file->celement=0;
  file->canode =0;
  file->cph_id =0;
  file->csrc_id=0;

  return(file);
}


void freeLADEventListFile(LADEventListFile** const file, 
			  int* const status)
{
  if (NULL!=*file) {
    if (NULL!=(*file)->fptr) {
      // If the file was opened in READWRITE mode, calculate
      // the check sum an append it to the FITS header.
      int mode;
      fits_file_mode((*file)->fptr, &mode, status);
      if (READWRITE==mode) {
	fits_write_chksum((*file)->fptr, status);
      }
      fits_close_file((*file)->fptr, status);
    }
    free(*file);
    *file=NULL;
  }
}


LADEventListFile* openNewLADEventListFile(const char* const filename,
					  const char clobber,
					  int* const status)
{
  LADEventListFile* file = newLADEventListFile(status);
  CHECK_STATUS_RET(*status, file);

  // Check if the file already exists.
  int exists;
  fits_file_exists(filename, &exists, status);
  CHECK_STATUS_RET(*status, file);
  if (0!=exists) {
    if (0!=clobber) {
      // Delete the file.
      remove(filename);
    } else {
      // Throw an error.
      char msg[MAXMSG];
      sprintf(msg, "file '%s' already exists", filename);
      SIXT_ERROR(msg);
      *status=EXIT_FAILURE;
      return(file);
    }
  }

  // Create a new LADEvent list FITS file from the template file.
  char buffer[MAXFILENAME];
  sprintf(buffer, "%s(%s%s)", filename, SIXT_DATA_PATH, 
	  "/templates/ladeventlist.tpl");
  fits_create_file(&file->fptr, buffer, status);
  CHECK_STATUS_RET(*status, file);

  // Set the time-keyword in the Event List Header.
  char datestr[MAXMSG];
  int timeref;
  fits_get_system_time(datestr, &timeref, status);
  CHECK_STATUS_RET(*status, file);
  fits_update_key(file->fptr, TSTRING, "DATE", datestr, 
		  "File creation date", status);
  CHECK_STATUS_RET(*status, file);

  // Add header information about program parameters.
  // The second parameter "1" means that the headers are written
  // to the first extension.
  HDpar_stamp(file->fptr, 1, status);
  CHECK_STATUS_RET(*status, file);

  // Move to the binary table extension.
  fits_movabs_hdu(file->fptr, 2, 0, status);
  CHECK_STATUS_RET(*status, file);

  // Close the file.
  freeLADEventListFile(&file, status);
  CHECK_STATUS_RET(*status, file);

  // Re-open the file.
  file = openLADEventListFile(filename, READWRITE, status);
  CHECK_STATUS_RET(*status, file);
  
  return(file);
}


LADEventListFile* openLADEventListFile(const char* const filename,
				       const int mode, 
				       int* const status)
{
  LADEventListFile* file = newLADEventListFile(status);
  CHECK_STATUS_RET(*status, file);

  headas_chat(4, "open event list file '%s' ...\n", filename);
  fits_open_table(&file->fptr, filename, mode, status);
  CHECK_STATUS_RET(*status, file);

  // Determine the row numbers.
  file->row=0;
  fits_get_num_rows(file->fptr, &file->nrows, status);
  CHECK_STATUS_RET(*status, file);

  // Determine the column numbers.
  fits_get_colnum(file->fptr, CASEINSEN, "TIME", &file->ctime, status);
  fits_get_colnum(file->fptr, CASEINSEN, "SIGNAL", &file->csignal, status);
  fits_get_colnum(file->fptr, CASEINSEN, "PANEL", &file->cpanel, status);
  fits_get_colnum(file->fptr, CASEINSEN, "MODULE", &file->cmodule, status);
  fits_get_colnum(file->fptr, CASEINSEN, "ELEMENT", &file->celement, status);
  fits_get_colnum(file->fptr, CASEINSEN, "ANODE", &file->canode, status);
  fits_get_colnum(file->fptr, CASEINSEN, "PH_ID", &file->cph_id, status);
  fits_get_colnum(file->fptr, CASEINSEN, "SRC_ID", &file->csrc_id, status);
  CHECK_STATUS_RET(*status, file);

  // Check if the vector length of the PH_ID and SRC_ID columns is equivalent 
  // with the corresponding array lengths in the LADEvent data structure.
  int typecode;
  long repeat, width;
  // PH_ID.
  fits_get_coltype(file->fptr, file->cph_id, &typecode, &repeat,
		   &width, status);
  CHECK_STATUS_RET(*status, file);
  if (repeat!=NLADEVENTPHOTONS) {
    // Throw an error.
    *status = EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "maximum number of photons contributing "
	    "to a single event is different "
	    "in the simulation (%d) and in the event list "
	    "template file (%ld)", NLADEVENTPHOTONS, repeat);
    SIXT_ERROR(msg);
    return(file);
  }
  // SRC_ID.
  fits_get_coltype(file->fptr, file->csrc_id, &typecode, &repeat,
		   &width, status);
  CHECK_STATUS_RET(*status, file);
  if (repeat!=NLADEVENTPHOTONS) {
    // Throw an error.
    *status = EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "maximum number of photons contributing "
	    "to a single event is different "
	    "in the simulation (%d) and in the event list "
	    "template file (%ld)!\n", NLADEVENTPHOTONS, repeat);
    SIXT_ERROR(msg);
    return(file);
  }

  return(file);
}


void addLADEvent2File(LADEventListFile* const file, 
		      LADEvent* const event, 
		      int* const status)
{
  // Check if the event file has been opened.
  CHECK_NULL_VOID(file, *status, "no event file opened");
  CHECK_NULL_VOID(file->fptr, *status, "no event file opened");

  // Insert a new, empty row to the table:
  fits_insert_rows(file->fptr, file->row, 1, status);
  CHECK_STATUS_VOID(*status);
  file->nrows++;
  file->row = file->nrows;

  // Write the data.
  updateLADEventInFile(file, file->row, event, status);
  CHECK_STATUS_VOID(*status);
}


void getLADEventFromFile(const LADEventListFile* const file,
			 const int row, LADEvent* const event,
			 int* const status)
{
  // Check if the file has been opened.
  CHECK_NULL_VOID(file, *status, "no event file opened");
  CHECK_NULL_VOID(file->fptr, *status, "no event file opened");

  // Check if there is still a row available.
  if (row > file->nrows) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("event file contains no further entries");
    return;
  }

  // Read in the data.
  int anynul=0;
  double dnull=0.;
  float fnull=0.;
  long lnull=0;

  fits_read_col(file->fptr, TDOUBLE, file->ctime, row, 1, 1, 
		&dnull, &event->time, &anynul, status);
  CHECK_STATUS_VOID(*status);

  fits_read_col(file->fptr, TFLOAT, file->csignal, row, 1, 1, 
		&fnull, &event->signal, &anynul, status);
  CHECK_STATUS_VOID(*status);

  fits_read_col(file->fptr, TLONG, file->cpanel, row, 1, 1, 
		&lnull, &event->panel, &anynul, status);
  CHECK_STATUS_VOID(*status);

  fits_read_col(file->fptr, TLONG, file->cmodule, row, 1, 1, 
		&lnull, &event->module, &anynul, status);
  CHECK_STATUS_VOID(*status);

  fits_read_col(file->fptr, TLONG, file->celement, row, 1, 1, 
		&lnull, &event->element, &anynul, status);
  CHECK_STATUS_VOID(*status);

  fits_read_col(file->fptr, TLONG, file->canode, row, 1, 1, 
		&lnull, &event->anode, &anynul, status);
  CHECK_STATUS_VOID(*status);

  fits_read_col(file->fptr, TLONG, file->cph_id, row, 1, NLADEVENTPHOTONS, 
		&lnull, &event->ph_id, &anynul, status);
  CHECK_STATUS_VOID(*status);

  fits_read_col(file->fptr, TLONG, file->csrc_id, row, 1, NLADEVENTPHOTONS, 
		&lnull, &event->src_id, &anynul, status);
  CHECK_STATUS_VOID(*status);

  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    *status = EXIT_FAILURE;
    SIXT_ERROR("reading from ImpactListFile failed");
    return;
  }
}


void updateLADEventInFile(const LADEventListFile* const file,
			  const int row, LADEvent* const event,
			  int* const status)
{
  fits_write_col(file->fptr, TDOUBLE, file->ctime, row, 
		 1, 1, &event->time, status);
  fits_write_col(file->fptr, TFLOAT, file->csignal, row, 
		 1, 1, &event->signal, status);
  fits_write_col(file->fptr, TLONG, file->cpanel, row, 
		 1, 1, &event->panel, status);
  fits_write_col(file->fptr, TLONG, file->cmodule, row, 
		 1, 1, &event->module, status);
  fits_write_col(file->fptr, TLONG, file->celement, row, 
		 1, 1, &event->element, status);
  fits_write_col(file->fptr, TLONG, file->canode, row, 
		 1, 1, &event->anode, status);
  fits_write_col(file->fptr, TLONG, file->cph_id, row, 
		 1, NLADEVENTPHOTONS, &event->ph_id, status);
  fits_write_col(file->fptr, TLONG, file->csrc_id, row, 
		 1, NLADEVENTPHOTONS, &event->src_id, status);
  CHECK_STATUS_VOID(*status);
}

