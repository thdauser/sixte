#include "ladimpactlistfile.h"


LADImpactListFile* newLADImpactListFile(int* const status)
{
  LADImpactListFile* file = (LADImpactListFile*)malloc(sizeof(LADImpactListFile));
  if (NULL==file) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for LADImpactListFile failed!\n", 
		   *status);
    return(file);
  }

  // Initialize pointers with NULL.
  file->fptr=NULL;

  // Initialize values.
  file->nrows=0;
  file->row  =0;
  file->ctime=0;
  file->cenergy =0;
  file->cpanel  =0;
  file->cmodule =0;
  file->celement=0;
  file->cx      =0;
  file->cy      =0;
  file->cph_id  =0;
  file->csrc_id =0;

  return(file);
}


void freeLADImpactListFile(LADImpactListFile** const file, int* const status)
{
  if (NULL!=*file) {
    if (NULL!=(*file)->fptr) {
      fits_close_file((*file)->fptr, status);
      headas_chat(5, "closed impact list file (containing %ld rows).\n", 
		  (*file)->nrows);
    }
    free(*file);
    *file=NULL;
  }
}


LADImpactListFile* openLADImpactListFile(const char* const filename,
					 const int mode, int* const status)
{
  LADImpactListFile* file = newLADImpactListFile(status);
  if (EXIT_SUCCESS!=*status) return(file);

  headas_chat(4, "open impact list file '%s' ...\n", filename);

  // Open the FITS file table for reading:
  if (fits_open_table(&file->fptr, filename, mode, status)) return(file);;

  // Get the HDU type.
  int hdutype;
  if (fits_get_hdu_type(file->fptr, &hdutype, status)) return(file);;
  // Image HDU results in an error message.
  if (IMAGE_HDU==hdutype) {
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "Error: no table extension available in FITS file '%s'!\n", 
	    filename);
    HD_ERROR_THROW(msg, *status);
    return(file);
  }

  // Determine the number of rows in the impact list.
  if (fits_get_num_rows(file->fptr, &file->nrows, status)) return(file);
  // Set internal row counter to the beginning of the file (starting at 0).
  file->row = 0;

  // Determine the individual column numbers.
  if(fits_get_colnum(file->fptr, CASEINSEN, "TIME", &file->ctime, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "ENERGY", &file->cenergy, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "PANEL", &file->cpanel, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "MODULE", &file->cmodule, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "ELEMENT", &file->celement, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "X", &file->cx, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "Y", &file->cy, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "PH_ID", &file->cph_id, status)) 
    return(file);
  if(fits_get_colnum(file->fptr, CASEINSEN, "SRC_ID", &file->csrc_id, status)) 
    return(file);

  return(file);
}


LADImpactListFile* openNewLADImpactListFile(const char* const filename,
					    int* const status)
{
  LADImpactListFile* file = newLADImpactListFile(status);
  if (EXIT_SUCCESS!=*status) return(file);

  // Remove old file, if it exists.
  remove(filename);

  // Create a new event list FITS file from the template file.
  char buffer[MAXFILENAME];
  sprintf(buffer, "%s(%s%s)", filename, SIXT_DATA_PATH, "/templates/ladimpactlist.tpl");
  fits_create_file(&file->fptr, buffer, status);
  CHECK_STATUS_RET(*status, file);

  // Set the time-keyword in the header.
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
	if (fits_update_key(file->fptr, TSTRING, "DATE-OBS", current_time_str, 
			    "Start Time (UTC) of exposure", status)) 
	  return(file);
      }
    }
  } 
  // END of writing time information to Event File FITS header.

  // Add header information about program parameters.
  // The second parameter "1" means that the headers are written
  // to the first extension.
  HDpar_stamp(file->fptr, 1, status);
  CHECK_STATUS_RET(*status, file);

  // Move to the binary table extension.
  fits_movabs_hdu(file->fptr, 2, 0, status);
  CHECK_STATUS_RET(*status, file);

  // Close the new LADImpactListFile.
  freeLADImpactListFile(&file, status);
  CHECK_STATUS_RET(*status, file);
  
  // Re-open the file.
  file = openLADImpactListFile(filename, READWRITE, status);
  CHECK_STATUS_RET(*status, file);

  return(file);
}


void getNextLADImpactFromFile(LADImpactListFile* const file, 
			      LADImpact* const impact, 
			      int* const status)
{
  // Check if the file has been opened.
  if (NULL==file) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: no LADImpactListFile opened!\n", *status);
    return;
  }
  if (NULL==file->fptr) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: no LADImpactListFile opened!\n", *status);
    return;
  }

  // Move counter to next line.
  file->row++;

  // Check if there is still a row available.
  if (file->row > file->nrows) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: LADImpactListFile contains no further entries!\n", *status);
    return;
  }

  // Read in the data.
  int anynul = 0;
  impact->time = 0.;
  fits_read_col(file->fptr, TDOUBLE, file->ctime, file->row, 1, 1, 
		&impact->time, &impact->time, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->energy = 0.;
  fits_read_col(file->fptr, TFLOAT, file->cenergy, file->row, 1, 1, 
		&impact->energy, &impact->energy, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->panel = 0;
  fits_read_col(file->fptr, TLONG, file->cpanel, file->row, 1, 1, 
		&impact->panel, &impact->panel, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->module = 0;
  fits_read_col(file->fptr, TLONG, file->cmodule, file->row, 1, 1, 
		&impact->module, &impact->module, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->element = 0;
  fits_read_col(file->fptr, TLONG, file->celement, file->row, 1, 1, 
		&impact->element, &impact->element, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->position.x = 0.;
  fits_read_col(file->fptr, TDOUBLE, file->cx, file->row, 1, 1, 
		&impact->position.x, &impact->position.x, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->position.y = 0.;
  fits_read_col(file->fptr, TDOUBLE, file->cy, file->row, 1, 1, 
		&impact->position.y, &impact->position.y, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->ph_id = 0;
  fits_read_col(file->fptr, TLONG, file->cph_id, file->row, 1, 1, 
		&impact->ph_id, &impact->ph_id, &anynul, status);
  CHECK_STATUS_VOID(*status);

  impact->src_id = 0;
  fits_read_col(file->fptr, TLONG, file->csrc_id, file->row, 1, 1, 
		&impact->src_id, &impact->src_id, &anynul, status);
  CHECK_STATUS_VOID(*status);
  
  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: reading from LADImpactListFile failed!\n", *status);
    return;
  }

  return;
}


void addLADImpact2File(LADImpactListFile* const ilf, 
		       LADImpact* const impact, 
		       int* const status)
{
  fits_insert_rows(ilf->fptr, ilf->row++, 1, status);
  fits_write_col(ilf->fptr, TDOUBLE, ilf->ctime, 
		 ilf->row, 1, 1, &impact->time, status);
  fits_write_col(ilf->fptr, TFLOAT, ilf->cenergy, 
		 ilf->row, 1, 1, &impact->energy, status);
  fits_write_col(ilf->fptr, TLONG, ilf->cpanel, 
		 ilf->row, 1, 1, &impact->panel, status);
  fits_write_col(ilf->fptr, TLONG, ilf->cmodule, 
		 ilf->row, 1, 1, &impact->module, status);
  fits_write_col(ilf->fptr, TLONG, ilf->celement, 
		 ilf->row, 1, 1, &impact->element, status);
  fits_write_col(ilf->fptr, TDOUBLE, ilf->cx, 
		 ilf->row, 1, 1, &(impact->position.x), status);
  fits_write_col(ilf->fptr, TDOUBLE, ilf->cy, 
		 ilf->row, 1, 1, &(impact->position.y), status);
  fits_write_col(ilf->fptr, TLONG, ilf->cph_id, 
		 ilf->row, 1, 1, &impact->ph_id, status);
  fits_write_col(ilf->fptr, TLONG, ilf->csrc_id, 
		 ilf->row, 1, 1, &impact->src_id, status);
  ilf->nrows++;
}

