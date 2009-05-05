#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif



#include "eventlist.h"



//////////////////////////////////////////////////////////////////
void add_eventlist_row(
		       struct Eventlist_File* eventlist_file, 
		       struct Event event, int *status
		       ) 
{
  fits_insert_rows(eventlist_file->fptr, eventlist_file->row++, 1, status);
  fits_write_col(eventlist_file->fptr, TDOUBLE, 1, eventlist_file->row, 
		 1, 1, &event.time, status);
  fits_write_col(eventlist_file->fptr, TLONG, 2, eventlist_file->row, 
		 1, 1, &event.pha, status);
  fits_write_col(eventlist_file->fptr, TINT, 3, eventlist_file->row, 
		 1, 1, &event.grade, status);
  fits_write_col(eventlist_file->fptr, TINT, 4, eventlist_file->row, 
		 1, 1, &event.xi, status);
  fits_write_col(eventlist_file->fptr, TINT, 5, eventlist_file->row, 
		 1, 1, &event.yi, status);
  fits_write_col(eventlist_file->fptr, TLONG, 6, eventlist_file->row, 
		 1, 1, &event.frame, status);

  // Set default values for PATNUM and PATID:
  // PATID has to be set to -1 !! 
  // Otherwise the pattern recognition algorithm doesn't work properly.
  // PATNUM
  event.patnum = 0;
  fits_write_col(eventlist_file->fptr, TLONG, 7, eventlist_file->row, 
		 1, 1, &event.patnum, status);  
  // PATID
  event.patid = -1;
  fits_write_col(eventlist_file->fptr, TLONG, 8, eventlist_file->row, 
		 1, 1, &event.patid, status);
  // Pile-up
  event.pileup = 0;
  fits_write_col(eventlist_file->fptr, TLONG, 9, eventlist_file->row, 
		 1, 1, &event.pileup, status); 
}



///////////////////////////////////////////////////////////////////
int create_eventlist_file(
			  struct Eventlist_File* eventlist_file,
			  Detector* detector,
			  double tstart,
			  double tend,
			  char *telescope_name,
			  char *ccd_name,
			  char *instrument_name,
			  int *status
			  )
{
  char *ftype[N_EVENT_FIELDS];
  char *fform[N_EVENT_FIELDS];
  char *funit[N_EVENT_FIELDS];
  int counter;

  char msg[MAXMSG];  // error output buffer

  do {   // Beginning of ERROR handling loop

    // Create a new FITS file:
    if (fits_create_file(&eventlist_file->fptr, eventlist_file->filename, status)) 
      break;

    // To create a FITS table, the format of the individual columns has to 
    // be specified.
    for(counter=0; counter<N_EVENT_FIELDS; counter++) {
      // Allocate memory
      ftype[counter] = (char *) malloc(8 * sizeof(char));
      fform[counter] = (char *) malloc(4 * sizeof(char));
      funit[counter] = (char *) malloc(10 * sizeof(char));
      
      // Check if all memory was allocated successfully:
      if ((!ftype[counter]) || (!fform[counter]) || (!funit[counter])) {
	*status = EXIT_FAILURE;
	sprintf(msg, "Error: no memory allocation for FITS table parameters "
		"failed (event list)!\n");
	HD_ERROR_THROW(msg, *status);
      }
    }

    // If an error has occurred during memory allocation, 
    // skip the following part.
    if (*status != EXIT_SUCCESS) break;

    // Set the field types of the table in the FITS file.
    // 1. time
    strcpy(ftype[0], "TIME");
    strcpy(fform[0], "D");
    strcpy(funit[0], "s");

    // 2. energy
    strcpy(ftype[1], "PHA");
    strcpy(fform[1], "J");
    strcpy(funit[1], "channel");

    // 3. grade
    strcpy(ftype[2], "GRADE");
    strcpy(fform[2], "I");
    strcpy(funit[2], "");

    // 4. x coordinate (integer pixel)
    strcpy(ftype[3], "COLUMN");   // RAWX
    strcpy(fform[3], "I");
    strcpy(funit[3], "");

    // 5. y coordinate (integer pixel)
    strcpy(ftype[4], "ROW");      // RAWY
    strcpy(fform[4], "I");
    strcpy(funit[4], "");

    // 6. frame number
    strcpy(ftype[5], "FRAME");
    strcpy(fform[5], "J");
    strcpy(funit[5], "");
  
    // 7. pattern number
    strcpy(ftype[6], "PATNUM");
    strcpy(fform[6], "J");
    strcpy(funit[6], "");

    // 8. pattern ID
    strcpy(ftype[7], "PATID");
    strcpy(fform[7], "J");
    strcpy(funit[7], "");

    // 9. pileup
    strcpy(ftype[8], "PILEUP");
    strcpy(fform[8], "J");
    strcpy(funit[8], "");


    // create the event list table
    if (fits_create_tbl(eventlist_file->fptr, BINARY_TBL, 0, N_EVENT_FIELDS, 
			ftype, fform, funit, "EVENTS", status)) break;
    

    // write descriptory data into the header of the FITS file
    if (fits_write_key(eventlist_file->fptr, TSTRING, "COMMENT", "EVENTLIST",
		       "content: eventlist of eROSITA simulation measurement", 
		       status)) break;

    // general mission headers
    if (fits_write_key(eventlist_file->fptr, TSTRING, "MISSION", "SpectrumXGamma", 
			"name of the mission", status)) break;
    if (fits_write_key(eventlist_file->fptr, TSTRING, "COMMENT", "DESCRIPT", 
			"eventlist file from the eROSITA telescope simulation",
			status)) break;

    // date and time headers
    char creation_date[30];
    int timeref;            // is 0, if returned time is in UTC
    if (fits_get_system_time(creation_date, &timeref, status)) break;
    if (fits_write_key(eventlist_file->fptr, TSTRING, "DATE", "2008-03-10",
		       "FITS file creation date (yyyy-mm-dd)", status)) break;
    if (fits_write_key(eventlist_file->fptr, TSTRING, "DATE-OBS", 
		       "2008-03-10T10:00:00",
		       "start time for the orbit", status)) break;
    if (fits_write_key(eventlist_file->fptr, TSTRING, "DATE-END", "",
		       "end time of the orbit", status)) break;
    double dbuffer = 0.;
    if (fits_write_key(eventlist_file->fptr, TDOUBLE, "MJDSTART", &dbuffer,
		       "start time of the orbit in Julian date format", 
		       status)) break;
    if (fits_write_key(eventlist_file->fptr, TDOUBLE, "MJDEND", &dbuffer,
		       "end time of the orbit in Julian date format", 
		       status)) break;
    dbuffer = 0.;
    if (fits_write_key(eventlist_file->fptr, TDOUBLE, "TIMEZERO", &dbuffer,
		       "Clock correction", status)) break;
    long lbuffer = 0;
    if (fits_write_key(eventlist_file->fptr, TLONG, "MJDREFI", &lbuffer,
		       "integer part of reference time", 
		       status)) break;
    if (fits_write_key(eventlist_file->fptr, TDOUBLE, "MJDREFF", &dbuffer,
		       "fractional part of reference time", status)) break;
    if (fits_write_key(eventlist_file->fptr, TDOUBLE, "TSTART", &tstart,
		       "start time of the orbit", status)) break;
    if (fits_write_key(eventlist_file->fptr, TDOUBLE, "TEND", &tend,
		       "start time of the orbit", 
		       status)) break;

    // Obligatory detector data
    if (fits_write_key (eventlist_file->fptr, TSTRING, "TELESCOP", telescope_name, 
			"name of the telescope", status)) break;
    if (fits_write_key (eventlist_file->fptr, TSTRING, "DETNAM", ccd_name, 
			"name of the detector", status)) break;
    if (fits_write_key (eventlist_file->fptr, TSTRING, "INSTRUME", instrument_name, 
			"name of the instrument", status)) break;

    // Determine the CCD mode (FRAMESTORE, DEPFET, ...)
    char data_mode[20];
    if (detector->type == FRAMESTORE) {
      strcpy(data_mode, "FRAMESTORE");
    } else if (detector->type == DEPFET) {
      strcpy(data_mode, "DEPFET");
    } else if (detector->type == TES) {
      strcpy(data_mode, "TES MCal");
    } else {strcpy(data_mode, "unknown");}
    if (fits_write_key (eventlist_file->fptr, TSTRING, "DATAMODE", data_mode, "", 
			status)) break;

    if (fits_write_key(eventlist_file->fptr, TSTRING, "FILTER", "none", "", status))
      break;
    if (fits_write_key(eventlist_file->fptr, TINT, "DETWIDTH", &detector->width, 
		       "width (number of pixels) of the detector", status)) 
      break;
    if (fits_write_key(eventlist_file->fptr, TINT, "DETCHANS", 
		       &detector->rmf->NumberChannels, 
		       "number of detector channels", status)) break;

    // instrument data
    if (fits_write_key (eventlist_file->fptr, TSTRING, "CREATOR", "simulation", 
			"", status)) break;

    // observation data
    if (fits_write_key (eventlist_file->fptr, TSTRING, "OBS_MODE", 
			"all-sky survey", "", status)) break;

    // Additional detector information:
    if (fits_write_key (eventlist_file->fptr, TSTRING, "COLUMN", ftype[3], 
			"column of event", status)) break;
    if (fits_write_key (eventlist_file->fptr, TSTRING, "ROW", ftype[4], 
			"row of event", status)) break;
    if (fits_write_key (eventlist_file->fptr, TINT, "COLUMNS", &detector->width, 
			"Number of columns", status)) break;
    if (fits_write_key (eventlist_file->fptr, TINT, "ROWS", &detector->width, 
			"Number of rows", status)) break;
    


    // If desired by the user, print all program parameters to HISTORY of 
    // FITS file (HDU number 1).
    HDpar_stamp(eventlist_file->fptr, 2, status);
    

    // Set initial value for new event list (start at beginning of the file).
    eventlist_file->row=0;  
    eventlist_file->nrows=0;

  } while (0);  // end of error handling loop


  //----------------
  // clean up
  for (counter=0; counter<N_EVENT_FIELDS; counter++) {
    if (ftype[counter]) free(ftype[counter]);
    if (fform[counter]) free(fform[counter]);
    if (funit[counter]) free(funit[counter]);
  }

  return (*status);
}




///////////////////////////////////////////////////////////////////////
// This function reads a row of data from the event list FITS table.
int get_eventlist_row(struct Eventlist_File eventlist_file, 
		     struct Event* event,
		     int *status) {
  int anynul = 0;

  // time (1st column)
  fits_read_col(eventlist_file.fptr, TDOUBLE, 1, eventlist_file.row+1, 1, 1, 
		&event->time, &event->time, &anynul, status);

  // PHA channel (2nd column)
  fits_read_col(eventlist_file.fptr, TLONG, 2, eventlist_file.row+1, 1, 1, 
		&event->pha, &event->pha, &anynul, status);

  // Grade (specifies split events)
  fits_read_col(eventlist_file.fptr, TINT, 3, eventlist_file.row+1, 1, 1, 
		&event->grade, &event->grade, &anynul, status);

  // xi (3rd column)
  fits_read_col(eventlist_file.fptr, TINT, 4, eventlist_file.row+1, 1, 1, 
		&event->xi, &event->xi, &anynul, status);

  // yi (4th column)
  fits_read_col(eventlist_file.fptr, TINT, 5, eventlist_file.row+1, 1, 1, 
		&event->yi, &event->yi, &anynul, status);

  // frame
  fits_read_col(eventlist_file.fptr, TLONG, 6, eventlist_file.row+1, 1, 1, 
		&event->frame, &event->frame, &anynul, status);

  // patnum
  fits_read_col(eventlist_file.fptr, TLONG, 7, eventlist_file.row+1, 1, 1, 
		&event->patnum, &event->patnum, &anynul, status);

  // patid
  fits_read_col(eventlist_file.fptr, TLONG, 8, eventlist_file.row+1, 1, 1, 
		&event->patid, &event->patid, &anynul, status);

  return(anynul);
}





///////////////////////////////////////////////////////////////////
// OBSOLETE !!
int open_eventlist_file(
			struct Eventlist_File* eventlist_file,
			int* status
			)
{
  char msg[MAXMSG];  // buffer for error messages

  printf("\nWarning: Using obsolete event list opening routine!\n");

  do {  // ERROR handling loop

    if (fits_open_table(&eventlist_file->fptr, eventlist_file->filename, 
			READONLY, status)) break;

    // get the HDU type
    int hdutype;
    if (fits_get_hdu_type(eventlist_file->fptr, &hdutype, status)) break;

    // image HDU results in an error message
    if (hdutype==IMAGE_HDU) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: no table extension available in event list "
	      "FITS file '%s'!\n", eventlist_file->filename);
      HD_ERROR_THROW(msg, *status);
      break;
    }

    // determine number of rows in the event list
    if (fits_get_num_rows(eventlist_file->fptr, &eventlist_file->nrows, status)) 
      break;

    // Set internal row counter to first row (starting at 0).
    eventlist_file->row = 0;

  } while(0);  // END of error handling loop

  return(*status);
}


///////////////////////////////////////////////////////////////////
// Opens an existing FITS file with a binary table event list
// for reading access.
int open_EventlistFile(char* filename, int* status)
{
  char msg[MAXMSG];  // buffer for error messages
  struct Eventlist_File *eventlistfile = NULL;
  
  do {  // ERROR handling loop
    
    // Memory allocation:
    eventlistfile = (struct Eventlist_File*)malloc(sizeof(struct Eventlist_File));
    if (NULL==eventlistfile) {
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: memory allocation in event list open routine failed!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }

    // Open the FITS file table for reading:
    if (fits_open_table(&eventlistfile->fptr, filename, 
			READONLY, status)) break;

    // Get the HDU type
    int hdutype;
    if (fits_get_hdu_type(eventlistfile->fptr, &hdutype, status)) break;

    // image HDU results in an error message
    if (hdutype==IMAGE_HDU) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: no table extension available in event list "
	      "FITS file '%s'!\n", filename);
      HD_ERROR_THROW(msg, *status);
      break;
    }

    // Determine the number of rows in the event list.
    if (fits_get_num_rows(eventlistfile->fptr, &eventlistfile->nrows, status)) 
      break;

    // Set internal row counter to first row (starting at 0).
    eventlistfile->row = 0;

  } while(0);  // END of error handling loop

  return(eventlistfile);
}

