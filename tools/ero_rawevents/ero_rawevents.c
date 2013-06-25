#include "ero_rawevents.h"


int ero_rawevents_main() 
{
  // Containing all programm parameters read by PIL
  struct Parameters par; 

  // Input event file.
  EventListFile* elf=NULL;

  // File pointer to the output eROSITA event file. 
  fitsfile* fptr=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL:
  set_toolname("ero_rawevents");
  set_toolversion("0.01");


  do { // Beginning of the ERROR handling loop.

    // --- Initialization ---

    headas_chat(3, "initialization ...\n");

    // Read parameters using PIL library:
    if ((status=getpar(&par))) break;


    // Open the input event file.
    elf=openEventListFile(par.EventList, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Read keywords from the input file.
    char comment[MAXMSG];
    float timezero=0.0;
    fits_read_key(elf->fptr, TFLOAT, "TIMEZERO", &timezero, comment, &status);
    CHECK_STATUS_BREAK(status);

    char date_obs[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "DATE-OBS", date_obs, comment, &status);
    CHECK_STATUS_BREAK(status);

    char time_obs[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "TIME-OBS", time_obs, comment, &status);
    CHECK_STATUS_BREAK(status);

    char date_end[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "DATE-END", date_end, comment, &status);
    CHECK_STATUS_BREAK(status);

    char time_end[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "TIME-END", time_end, comment, &status);
    CHECK_STATUS_BREAK(status);

    double tstart=0.0;
    fits_read_key(elf->fptr, TDOUBLE, "TSTART", &tstart, comment, &status);
    CHECK_STATUS_BREAK(status);

    double tstop=0.0;
    fits_read_key(elf->fptr, TDOUBLE, "TSTOP", &tstop, comment, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the file creation date for the header.
    char creation_date[MAXMSG];
    int timeref;
    fits_get_system_time(creation_date, &timeref, &status);
    CHECK_STATUS_BREAK(status);

    // Check if the output file already exists.
    int exists;
    fits_file_exists(par.eroEventList, &exists, &status);
    CHECK_STATUS_BREAK(status);
    if (0!=exists) {
      if (0!=par.clobber) {
	// Delete the file.
	remove(par.eroEventList);
      } else {
	// Throw an error.
	char msg[MAXMSG];
	sprintf(msg, "file '%s' already exists", par.eroEventList);
	SIXT_ERROR(msg);
	status=EXIT_FAILURE;
	break;
      }
    }
    
    // Create and open a new FITS file.
    headas_chat(3, "create new eROSITA event list file '%s' ...\n",
		par.eroEventList);
    fits_create_file(&fptr, par.eroEventList, &status);
    CHECK_STATUS_BREAK(status);

    // Create the event table.
    char *ttype[]={"TIME", "FRAME", "RAWX", "RAWY", "PHA"};
    char *tunit[]={"", "", "", "", "ADU"};
    char *tform[]={"D", "J", "I", "I", "I"}; 
    fits_create_tbl(fptr, BINARY_TBL, 0, 5, ttype, tform, tunit, 
		    "EVENTS", &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not create binary table for events "
	      "in file '%s'", par.eroEventList);
      SIXT_ERROR(msg);
      break;
    }

    // Insert header keywords.
    char hduclass[MAXMSG]="OGIP";
    fits_update_key(fptr, TSTRING, "HDUCLASS", hduclass, "", &status);
    char hduclas1[MAXMSG]="EVENTS";
    fits_update_key(fptr, TSTRING, "HDUCLAS1", hduclas1, "", &status);
    CHECK_STATUS_BREAK(status);

    // Insert the standard eROSITA header keywords.
    sixt_add_fits_erostdkeywords(fptr, 1, creation_date, date_obs, time_obs,
				 date_end, time_end, tstart, tstop, 
				 timezero, &status);
    CHECK_STATUS_BREAK(status);
    sixt_add_fits_erostdkeywords(fptr, 2, creation_date, date_obs, time_obs,
				 date_end, time_end, tstart, tstop, 
				 timezero, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the column numbers.
    int ctime, cframe, crawx, crawy, cpha;
    fits_get_colnum(fptr, CASEINSEN, "TIME", &ctime, &status);
    fits_get_colnum(fptr, CASEINSEN, "FRAME", &cframe, &status);
    fits_get_colnum(fptr, CASEINSEN, "RAWX", &crawx, &status);
    fits_get_colnum(fptr, CASEINSEN, "RAWY", &crawy, &status);
    fits_get_colnum(fptr, CASEINSEN, "PHA", &cpha, &status);
    CHECK_STATUS_BREAK(status);

    // Set the TLMIN and TLMAX keywords for the PHA column.
    char keyword[MAXMSG];
    int tlmin_pha=0, tlmax_pha=4095;
    sprintf(keyword, "TLMIN%d", cpha);
    fits_update_key(fptr, TINT, keyword, &tlmin_pha, "", &status);
    sprintf(keyword, "TLMAX%d", cpha);
    fits_update_key(fptr, TINT, keyword, &tlmax_pha, "", &status);
    CHECK_STATUS_BREAK(status);

    // --- END of Initialization ---

    
    // --- Beginning of Copy Process ---

    headas_chat(3, "start copy process ...\n");

    // Loop over all events in the input file. 
    long row;
    for (row=0; row<elf->nrows; row++) {
      
      // Read the next event from the input file.
      Event event;
      getEventFromFile(elf, row+1, &event, &status);
      CHECK_STATUS_BREAK(status);

      // Store it in the output file.
      fits_insert_rows(fptr, row, 1, &status);
      CHECK_STATUS_BREAK(status);

      fits_write_col(fptr, TDOUBLE, ctime, row+1, 1, 1, 
		     &event.time, &status);
      fits_write_col(fptr, TLONG, cframe, row+1, 1, 1, 
		     &event.frame, &status);
      fits_write_col(fptr, TLONG, cpha, row+1, 1, 1, 
		     &event.pi, &status);
      int rawx=event.rawx+1;
      fits_write_col(fptr, TINT, crawx, row+1, 1, 1, &rawx, &status);
      int rawy=event.rawy+1;
      fits_write_col(fptr, TINT, crawy, row+1, 1, 1, &rawy, &status);

      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);
    // END of loop over all events in the input file.

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Close the files.
  freeEventListFile(&elf, &status);

  // Append the check sum to the FITS header.
  if (NULL!=fptr) {
    int ii;
    for (ii=0; ii<2; ii++) {
      int hdutype;
      fits_movabs_hdu(fptr, ii+1, &hdutype, &status);
      fits_write_chksum(fptr, &status);
    }
    fits_close_file(fptr, &status);
  }
  
  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully\n\n");
  return(status);
}


int getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  status=ape_trad_query_file_name("EventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the input file");
    return(status);
  } 
  strcpy(par->EventList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("eroEventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the output file");
    return(status);
  } 
  strcpy(par->eroEventList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  return(status);
}

