#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include "eventlist.h"


//////////////////////////////////////////////////////////////////
void add_eventlist_row(struct Eventlist_File* ef, struct Event event, int *status) 
{
  // Insert a new, empty row to the table:
  fits_insert_rows(ef->fptr, ef->row++, 1, status);

  if (0<ef->ctime) 
    fits_write_col(ef->fptr, TDOUBLE, ef->ctime, ef->row, 1, 1, &event.time, status);
  if (0<ef->cpha)
    fits_write_col(ef->fptr, TLONG, ef->cpha, ef->row, 1, 1, &event.pha, status);
  if (0<ef->crawx) {
    event.xi += ef->PixelOffset;
    fits_write_col(ef->fptr, TINT, ef->crawx, ef->row, 1, 1, &event.xi, status);
  }
  if (0<ef->crawy) {
    event.yi += ef->PixelOffset;
    fits_write_col(ef->fptr, TINT, ef->crawy, ef->row, 1, 1, &event.yi, status);
  }
  if (0<ef->cframe)
    fits_write_col(ef->fptr, TLONG, ef->cframe, ef->row, 1, 1, &event.frame, status);

  // Set default values for PATNUM and PATID:
  // PATID has to be set to -1 !! 
  // Otherwise the pattern recognition algorithm doesn't work properly.

  // PATNUM
  if (0<ef->cpatnum) {
    event.patnum = 0;
    fits_write_col(ef->fptr, TLONG, ef->cpatnum, ef->row, 1, 1, &event.patnum, status);  
  }

  // PATID
  if (0<ef->cpatid) {
    event.patid = -1;
    fits_write_col(ef->fptr, TLONG, ef->cpatid, ef->row, 1, 1, &event.patid, status);
  }

  // Pile-up
  if (0<ef->cpileup) {
    event.pileup = 0;
    fits_write_col(ef->fptr, TLONG, ef->cpileup, ef->row, 1, 1, &event.pileup, status); 
  }

  // RA and DEC
  if (0<ef->cra) 
    fits_write_col(ef->fptr, TDOUBLE, ef->cra, ef->row, 1, 1, &event.ra, status);
  if (0<ef->cdec) 
    fits_write_col(ef->fptr, TDOUBLE, ef->cdec, ef->row, 1, 1, &event.dec, status);

  // Sky coordinates X and Y
  if (0<ef->cskyx) 
    fits_write_col(ef->fptr, TLONG, ef->cskyx, ef->row, 1, 1, &event.sky_xi, status);
  if (0<ef->cskyy) 
    fits_write_col(ef->fptr, TLONG, ef->cskyy, ef->row, 1, 1, &event.sky_yi, status);

}


/*
///////////////////////////////////////////////////////////////////
struct Eventlist_File* create_Eventlist_File(
					     char* filename,
					     Detector* detector,
					     double tstart,
					     double tend,
					     int *status
					     )
{
  struct Eventlist_File* ef=NULL;

  char mission[MAXMSG];
  char telescope[MAXMSG];
  char instrument[MAXMSG];
  char detname[MAXMSG];
  char data_mode[MAXMSG];

  char msg[MAXMSG];  // error output buffer

  do {   // Beginning of ERROR handling loop
    
    // Allocate memory
    ef = (struct Eventlist_File*)malloc(sizeof(struct Eventlist_File));
    if (NULL==ef) {
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: memory allocation for Eventlist_File data structure failed!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }
    struct Eventlist_File empty_list = { .fptr=NULL };
    *ef = empty_list;

    // Create a new FITS file:
    if (fits_create_file(&ef->fptr, filename, status)) break;

    // To create a FITS table, the format of the individual columns has to 
    // be specified.
    // Distingusih between the different detectors and determine the required 
    // number of table columns.
    switch (detector->type) {
    case FRAMESTORE: 
      ef->ncolumns = 9;
      break;
    case DEPFET:
      ef->ncolumns = 8;
      break;
    default:
      ef->ncolumns = 0;
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: invalid detector type in event list creation!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }
    if (EXIT_SUCCESS!=*status) break;
    ef->detectortype = detector->type;

    char *ftype[ef->ncolumns];
    char *fform[ef->ncolumns];
    char *funit[ef->ncolumns];
    
    // Allocate memory    
    int counter;
    for(counter=0; counter<ef->ncolumns; counter++) {
      ftype[counter] = (char *) malloc(10 * sizeof(char));
      fform[counter] = (char *) malloc(10 * sizeof(char));
      funit[counter] = (char *) malloc(10 * sizeof(char));
      
      // Check if all memory was allocated successfully:
      if ((!ftype[counter]) || (!fform[counter]) || (!funit[counter])) {
	*status = EXIT_FAILURE;
	sprintf(msg, "Error: no memory allocation for FITS table parameters "
		"failed (event list)!\n");
	HD_ERROR_THROW(msg, *status);
	break;
      }
    }
    // If an error has occurred during memory allocation, 
    // skip the following part.
    if (*status != EXIT_SUCCESS) break;


    // Define the Instrument-specific layout of the event list table:
    int column_counter = 0;
    switch (ef->detectortype) {
    case FRAMESTORE: // eROSITA pnCCD
      // Readout time
      strcpy(ftype[column_counter], "TIME"); 
      strcpy(fform[column_counter], "D"); 
      strcpy(funit[column_counter], "s");
      ef->ctime = ++column_counter;
      // PHA channel
      strcpy(ftype[column_counter], "PHA"); 
      strcpy(fform[column_counter], "J"); 
      strcpy(funit[column_counter], "channel");
      ef->cpha = ++column_counter;
      // RAWX and RAWY
      strcpy(ftype[column_counter], "RAWX"); 
      strcpy(fform[column_counter], "I"); 
      strcpy(funit[column_counter], "");
      ef->crawx = ++column_counter;
      strcpy(ftype[column_counter], "RAWY"); 
      strcpy(fform[column_counter], "I"); 
      strcpy(funit[column_counter], "");
      ef->crawy = ++column_counter;
      // FRAME
      strcpy(ftype[column_counter], "FRAME"); 
      strcpy(fform[column_counter], "J"); 
      strcpy(funit[column_counter], "");
      ef->cframe = ++column_counter;
      // Source position in equatorial coordinates RA and DEC
      strcpy(ftype[column_counter], "RA"); 
      strcpy(fform[column_counter], "D"); // R*8
      strcpy(funit[column_counter], "degree");
      ef->cra = ++column_counter;
      strcpy(ftype[column_counter], "DEC"); 
      strcpy(fform[column_counter], "D"); // R*8
      strcpy(funit[column_counter], "degree");
      ef->cdec = ++column_counter;
      // Sky coordinates in integer pixels of size 0.05"
      strcpy(ftype[column_counter], "X"); 
      strcpy(fform[column_counter], "J"); // I*4
      strcpy(funit[column_counter], "pixel");
      ef->cskyx = ++column_counter;
      strcpy(ftype[column_counter], "Y"); 
      strcpy(fform[column_counter], "J"); // I*4
      strcpy(funit[column_counter], "pixel");
      ef->cskyy = ++column_counter;

      ef->cpatid=0; ef->cpatnum=0; ef->cpileup=0;

      break;

    case DEPFET: // IXO WFI
      // Readout time
      strcpy(ftype[column_counter], "TIME"); 
      strcpy(fform[column_counter], "D"); 
      strcpy(funit[column_counter], "s");
      ef->ctime = ++column_counter;
      // PHA channel
      strcpy(ftype[column_counter], "PHA"); 
      strcpy(fform[column_counter], "J"); 
      strcpy(funit[column_counter], "channel");
      ef->cpha = ++column_counter;
      // RAWX and RAWY
      strcpy(ftype[column_counter], "COLUMN"); 
      strcpy(fform[column_counter], "I"); 
      strcpy(funit[column_counter], "");
      ef->crawx = ++column_counter;
      strcpy(ftype[column_counter], "ROW"); 
      strcpy(fform[column_counter], "I"); 
      strcpy(funit[column_counter], "");
      ef->crawy = ++column_counter;
      // FRAME
      strcpy(ftype[column_counter], "FRAME"); 
      strcpy(fform[column_counter], "J"); 
      strcpy(funit[column_counter], "");
      ef->cframe = ++column_counter;
      // Pattern number
      strcpy(ftype[column_counter], "PATNUM"); 
      strcpy(fform[column_counter], "J"); 
      strcpy(funit[column_counter], "");
      ef->cpatnum = ++column_counter;
      // Pattern ID
      strcpy(ftype[column_counter], "PATID"); 
      strcpy(fform[column_counter], "J"); 
      strcpy(funit[column_counter], "");
      ef->cpatid = ++column_counter;
      // Pileup
      strcpy(ftype[column_counter], "PILEUP"); 
      strcpy(fform[column_counter], "J"); 
      strcpy(funit[column_counter], "");
      ef->cpileup = ++column_counter;

      ef->cra=0; ef->cdec=0; ef->cskyx=0; ef->cskyy=0;

      break;

    default:
      break;
    }

    // Create the event list table in the FITS file.
    if (fits_create_tbl(ef->fptr, BINARY_TBL, 0, ef->ncolumns, 
			ftype, fform, funit, "EVENTS", status)) break;
    

    // HEADER Keywords:
    // Distinguish between the different missions:
    switch (ef->detectortype) {
    case FRAMESTORE:    // eROSITA
      strcpy(mission,   "Simulation"); // SRG
      strcpy(telescope, "FM4");
      strcpy(detname,   "pnCCD1" );
      strcpy(instrument,"eROSITA"); 
      strcpy(data_mode, "FRAMESTORE");
      break;

    case DEPFET:        // IXO WFI
      strcpy(mission,   "IXO");
      strcpy(telescope, "IXO");
      strcpy(detname,   "WFI" );
      strcpy(instrument,"WFI"); 
      strcpy(data_mode, "DEPFET");
      break;

    default:
      break;
    }

    // Write descriptory data into the HEADER of the FITS file.
    if (fits_write_key(ef->fptr, TSTRING, "COMMENT", "EVENTLIST", "", status)) break;

    // General mission headers
    if (fits_write_key(ef->fptr, TSTRING, "MISSION",  mission, 
		       "name of the mission", status)) break;
    if (fits_write_key(ef->fptr, TSTRING, "COMMENT", "DESCRIPT", 
		       "content: event list generated by SIXT "
		       "(SImulation of X-ray Telescopes)", status)) break;

    // Obligatory detector data
    if (fits_write_key (ef->fptr, TSTRING, "TELESCOP", telescope, 
			"name of the telescope", status)) break;
    if (fits_write_key (ef->fptr, TSTRING, "DETNAM", detname, 
			"name of the detector", status)) break;
    if (fits_write_key (ef->fptr, TSTRING, "INSTRUME", instrument, 
			"name of the instrument", status)) break;

    // Date and time headers
    char creation_date[30];
    int timeref;            // is 0, if returned time is in UTC
    if (fits_get_system_time(creation_date, &timeref, status)) break;
    if (fits_write_key(ef->fptr, TSTRING, "DATE", "2008-03-10",
		       "FITS file creation date (yyyy-mm-dd)", status)) break;
    if (fits_write_key(ef->fptr, TSTRING, "DATE-OBS", "2008-03-10T10:00:00",
		       "start time for the orbit", status)) break;
    if (fits_write_key(ef->fptr, TSTRING, "DATE-END", "",
		       "end time of the orbit", status)) break;
    double dbuffer = 0.;
    if (fits_write_key(ef->fptr, TDOUBLE, "MJDSTART", &dbuffer,
		       "start time of the orbit in Julian date format", 
		       status)) break;
    if (fits_write_key(ef->fptr, TDOUBLE, "MJDEND", &dbuffer,
		       "end time of the orbit in Julian date format", 
		       status)) break;
    dbuffer = 0.;
    if (fits_write_key(ef->fptr, TDOUBLE, "TIMEZERO", &dbuffer,
		       "Clock correction", status)) break;
    long lbuffer = 0;
    if (fits_write_key(ef->fptr, TLONG, "MJDREFI", &lbuffer,
		       "integer part of reference time", status)) break;
    if (fits_write_key(ef->fptr, TDOUBLE, "MJDREFF", &dbuffer,
		       "fractional part of reference time", status)) break;
    if (fits_write_key(ef->fptr, TDOUBLE, "TSTART", &tstart,
		       "start time of the orbit", status)) break;
    if (fits_write_key(ef->fptr, TDOUBLE, "TEND", &tend, 
		       "start time of the orbit", status)) break;


    // Determine the CCD mode (FRAMESTORE, DEPFET, ...)
    if (fits_write_key (ef->fptr, TSTRING, "DATAMODE", data_mode, "", status)) break;
    if (fits_write_key(ef->fptr, TSTRING, "FILTER", "none", "", status)) break;
    if (fits_write_key(ef->fptr, TINT, "DETWIDTH", &detector->width, 
		       "width (number of pixels) of the detector", status)) break;
    if (fits_write_key(ef->fptr, TINT, "DETCHANS", &detector->rmf->NumberChannels,
		       "number of detector channels", status)) break;

    // Instrument data
    if (fits_write_key (ef->fptr, TSTRING, "CREATOR", "SIXT", "", status)) break;

    // Observation data
    if (fits_write_key (ef->fptr, TSTRING, "OBS_MODE", "all-sky survey", "", status)) break;

    // Additional detector information ONLY for IXO WFI:
    if (DEPFET == ef->detectortype) {
      if (fits_write_key (ef->fptr, TSTRING, "COLUMN", ftype[3],
			  "column of event", status)) break;
      if (fits_write_key (ef->fptr, TSTRING, "ROW", ftype[4],
			  "row of event", status)) break;
      if (fits_write_key (ef->fptr, TINT, "COLUMNS", &detector->width, 
			  "Number of columns", status)) break;
      if (fits_write_key (ef->fptr, TINT, "ROWS", &detector->width, 
			  "Number of rows", status)) break;
    }
    

    // If desired by the user, print all program parameters to HISTORY of 
    // FITS file (HDU number 1).
    HDpar_stamp(ef->fptr, 2, status);
    

    // Set initial value for new event list (start at beginning of the file).
    ef->row=0;  
    ef->nrows=0;

    // clean up
    for (counter=0; counter<ef->ncolumns; counter++) {
      if (ftype[counter]) free(ftype[counter]);
      if (fform[counter]) free(fform[counter]);
      if (funit[counter]) free(funit[counter]);
    }

  } while (0);  // end of error handling loop

  return(ef);
}
*/







///////////////////////////////////////////////////////////////////////
// This function reads a row of data from the event list FITS table.
int get_eventlist_row(struct Eventlist_File ef, 
		      struct Event* event, int *status) 
{
  int anynul = 0;

  // time (1st column)
  event->time = 0.;
  if (0<ef.ctime)
    fits_read_col(ef.fptr, TDOUBLE, ef.ctime, ef.row+1, 1, 1, 
		  &event->time, &event->time, &anynul, status);

  // PHA channel (2nd column)
  event->pha = 0;
  if (0<ef.cpha)
    fits_read_col(ef.fptr, TLONG, ef.cpha, ef.row+1, 1, 1, 
		  &event->pha, &event->pha, &anynul, status);

  // xi (3rd column) (RAWX)
  event->xi = 0;
  if (0<ef.crawx)
    fits_read_col(ef.fptr, TINT, ef.crawx, ef.row+1, 1, 1, 
		  &event->xi, &event->xi, &anynul, status);
  event->xi -= ef.PixelOffset;

  // yi (4th column) (RAWY)
  event->yi = 0;
  if (0<ef.crawy)
    fits_read_col(ef.fptr, TINT, ef.crawy, ef.row+1, 1, 1, 
		  &event->yi, &event->yi, &anynul, status);
  event->yi -= ef.PixelOffset;

  // frame
  event->frame = 0;
  if (0<ef.cframe) 
    fits_read_col(ef.fptr, TLONG, ef.cframe, ef.row+1, 1, 1, 
		  &event->frame, &event->frame, &anynul, status);

  // patnum
  event->patnum = 0;
  if (0<ef.cpatnum)
    fits_read_col(ef.fptr, TLONG, ef.cpatnum, ef.row+1, 1, 1, 
		  &event->patnum, &event->patnum, &anynul, status);

  // patid
  event->patid = 0;
  if (0<ef.cpatid)
    fits_read_col(ef.fptr, TLONG, ef.cpatid, ef.row+1, 1, 1, 
		  &event->patid, &event->patid, &anynul, status);

  return(anynul);
}






///////////////////////////////////////////////////////////////////
struct Eventlist_File* open_EventlistFile(char* filename, int access_mode, int* status)
{
  char msg[MAXMSG];  // buffer for error messages
  struct Eventlist_File *ef = NULL;
  
  do {  // ERROR handling loop
    
    // Memory allocation:
    ef = (struct Eventlist_File*)malloc(sizeof(struct Eventlist_File));
    if (NULL==ef) {
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: memory allocation in event list open routine failed!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }

    // Open the FITS file table for reading:
    if (fits_open_table(&ef->fptr, filename, access_mode, status)) break;

    // Get the HDU type
    int hdutype;
    if (fits_get_hdu_type(ef->fptr, &hdutype, status)) break;

    // Image HDU results in an error message.
    if (IMAGE_HDU==hdutype) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: no table extension available in event list "
	      "FITS file '%s'!\n", filename);
      HD_ERROR_THROW(msg, *status);
      break;
    }

    // Determine the number of rows in the event list.
    if (fits_get_num_rows(ef->fptr, &ef->nrows, status)) break;

    // Set internal row counter to first row (starting at 0).
    ef->row = 0;


    // Determine the individual column numbers:
    // REQUIRED columns:
    int status2 = EXIT_SUCCESS;
    if(fits_get_colnum(ef->fptr, CASEINSEN, "TIME", &ef->ctime, status)) break;
    if(fits_get_colnum(ef->fptr, CASEINSEN, "PHA", &ef->cpha, status)) break;
    if(fits_get_colnum(ef->fptr, CASEINSEN, "FRAME", &ef->cframe, status)) break;

    if((fits_get_colnum(ef->fptr, CASEINSEN, "RAWX", &ef->crawx, status)) && 
       (fits_get_colnum(ef->fptr, CASEINSEN, "COLUMN", &ef->crawx, &status2))) {
      *status += status2; 
      break;
    } else {
      *status = EXIT_SUCCESS;
      status2 = EXIT_SUCCESS;
    }
    if((fits_get_colnum(ef->fptr, CASEINSEN, "RAWY", &ef->crawy, status)) && 
       (fits_get_colnum(ef->fptr, CASEINSEN, "ROW", &ef->crawy, &status2))) {
      *status += status2; 
      break;
    } else {
      *status = EXIT_SUCCESS;
      status2 = EXIT_SUCCESS;
    }


    // OPTIONAL columns:
    int opt_status = EXIT_SUCCESS;
    // eROSITA:
    opt_status = EXIT_SUCCESS;
    if(fits_get_colnum(ef->fptr, CASEINSEN, "RA", &ef->cra, &opt_status)) ef->cra=0;
    opt_status=EXIT_SUCCESS;
    if(fits_get_colnum(ef->fptr, CASEINSEN, "DEC", &ef->cdec, &opt_status)) ef->cdec=0;
    opt_status=EXIT_SUCCESS;
    if(fits_get_colnum(ef->fptr, CASEINSEN, "X", &ef->cskyx, &opt_status)) ef->cskyx=0;
    opt_status=EXIT_SUCCESS;
    if(fits_get_colnum(ef->fptr, CASEINSEN, "Y", &ef->cskyy, &opt_status)) ef->cskyy=0;

    // IXO WFI:
    opt_status=0;
    if(fits_get_colnum(ef->fptr, CASEINSEN, "PATNUM", &ef->cpatnum, &opt_status)) 
      ef->cpatnum=0;
    opt_status=EXIT_SUCCESS;
    if(fits_get_colnum(ef->fptr, CASEINSEN, "PATID", &ef->cpatid, &opt_status)) 
      ef->cpatid=0;
    opt_status=EXIT_SUCCESS;
    if(fits_get_colnum(ef->fptr, CASEINSEN, "PILEUP", &ef->cpileup, &opt_status)) 
      ef->cpileup=0;

    // Determine the PixelOffset (numbering scheme of RAWX and RAWY).
    // Mainly used for eROSITA, as the numbering starts at 1 instead of 0 there.
    // Default value is 0.
    char comment[MAXMSG];
    opt_status=EXIT_SUCCESS;
    if(fits_read_key(ef->fptr, TINT, "PXOFFSET", &ef->PixelOffset, comment, &opt_status)) {
      ef->PixelOffset = 0;
    }

    // Clear the HEAdas error stack in order to delete error messages created
    // due to missing OPTIONAL column names.
    HDerror_reset();

  } while(0);  // END of error handling loop

  if(EXIT_SUCCESS!=*status) ef=NULL;
  return(ef);
}

