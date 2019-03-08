/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "ero_rawevents.h"


int ero_rawevents_main()
{
  // Containing all programm parameters read by PIL
  struct Parameters par;

  // Input event file.
  EventFile* elf=NULL;

  // File pointer to the output eROSITA event file.
  fitsfile* fptr=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL:
  set_toolname("ero_rawevents");
  set_toolversion("0.03");


  do { // Beginning of the ERROR handling loop.

    // --- Initialization ---

    headas_chat(3, "initialization ...\n");

    // Read parameters using PIL library:
    if ((status=getpar(&par))) break;


    // Open the input event file.
    elf=openEventFile(par.RawData, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Check if the input file contains single-pixel events.
    char evtype[MAXMSG], comment[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "EVTYPE", evtype, comment, &status);
    if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("could not read FITS keyword 'EVTYPE'");
      break;
    }
    strtoupper(evtype);
    if (0!=strcmp(evtype, "PIXEL")) {
      status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "event type of input file is '%s' (must be 'PIXEL')", evtype);
      SIXT_ERROR(msg);
      break;
    }

    // Read keywords from the input file.
    double mjdref=0.0;
    fits_read_key(elf->fptr, TDOUBLE, "MJDREF", &mjdref, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'MJDREF' from input "
	      "event list '%s'", par.RawData);
      SIXT_ERROR(msg);
      break;
    }

    double timezero=0.0;
    fits_write_errmark();
    int status2=EXIT_SUCCESS;
    fits_read_key(elf->fptr, TDOUBLE, "TIMEZERO", &timezero, comment, &status2);
    fits_clear_errmark();
    if (EXIT_SUCCESS!=status2) {
      timezero=0.;
    }

    char date_obs[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "DATE-OBS", date_obs, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'DATE-OBS' from input "
	      "event list '%s'", par.RawData);
      SIXT_ERROR(msg);
      break;
    }

    char time_obs[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "TIME-OBS", time_obs, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'TIME-OBS' from input "
	      "event list '%s'", par.RawData);
      SIXT_ERROR(msg);
      break;
    }

    char date_end[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "DATE-END", date_end, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'DATE-END' from input "
	      "event list '%s'", par.RawData);
      SIXT_ERROR(msg);
      break;
    }

    char time_end[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "TIME-END", time_end, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'TIME-END' from input "
	      "event list '%s'", par.RawData);
      SIXT_ERROR(msg);
      break;
    }

    double tstart=0.0;
    fits_read_key(elf->fptr, TDOUBLE, "TSTART", &tstart, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'TSTART' from input "
	      "event list '%s'", par.RawData);
      SIXT_ERROR(msg);
      break;
    }

    double tstop=0.0;
    fits_read_key(elf->fptr, TDOUBLE, "TSTOP", &tstop, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'TSTOP' from input "
	      "event list '%s'", par.RawData);
      SIXT_ERROR(msg);
      break;
    }

    // Read the FILTER key word (from event file)
    char filter[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "FILTER", filter, comment, &status);
    if (EXIT_SUCCESS!=status) {
    	char msg[MAXMSG];
    	sprintf(msg, "could not read FITS keyword 'FILTER' from input "
    			"event list '%s'", par.RawData);
    	SIXT_ERROR(msg);
    	break;
    }

    // Verify values of MJDREF and TIMEZERO.
    verifyMJDREF(eromjdref, mjdref, "in event file", &status);
    CHECK_STATUS_BREAK(status);
    verifyTIMEZERO(timezero, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the file creation date for the header.
    char creation_date[MAXMSG];
    int timeref;
    fits_get_system_time(creation_date, &timeref, &status);
    CHECK_STATUS_BREAK(status);

    // Check if the output file already exists.
    int exists;
    fits_file_exists(par.eroEvtFile, &exists, &status);
    CHECK_STATUS_BREAK(status);
    if (0!=exists) {
      if (0!=par.clobber) {
	// Delete the file.
	remove(par.eroEvtFile);
      } else {
	// Throw an error.
	char msg[MAXMSG];
	sprintf(msg, "file '%s' already exists", par.eroEvtFile);
	SIXT_ERROR(msg);
	status=EXIT_FAILURE;
	break;
      }
    }

    // Create and open a new FITS file.
    headas_chat(3, "create new eROSITA event list file '%s' ...\n",
		par.eroEvtFile);
    fits_create_file(&fptr, par.eroEvtFile, &status);
    CHECK_STATUS_BREAK(status);

    // Create the event table.
    char *ttype[]={"OTS", "FRACSEC", "FRAME", "RAWX", "RAWY", "PHA"};
    char *tunit[]={"s", "micros", "", "", "", "ADU"};
    char *tform[]={"J", "J", "J", "I", "I", "I"};
    fits_create_tbl(fptr, BINARY_TBL, 0, 6, ttype, tform, tunit,
		    "EVENTS", &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not create binary table for events "
	      "in file '%s'", par.eroEvtFile);
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
    sixt_add_fits_erostdkeywords(fptr, 1, filter, creation_date, date_obs, time_obs,
				 date_end, time_end, tstart, tstop,
				 mjdref, timezero, par.CCDNr, &status);
    CHECK_STATUS_BREAK(status);
    sixt_add_fits_erostdkeywords(fptr, 2, filter, creation_date, date_obs, time_obs,
				 date_end, time_end, tstart, tstop,
				 mjdref, timezero, par.CCDNr, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the column numbers.
    int cots, cfracsec, cframe, crawx, crawy, cpha;
    fits_get_colnum(fptr, CASEINSEN, "OTS", &cots, &status);
    fits_get_colnum(fptr, CASEINSEN, "FRACSEC", &cfracsec, &status);
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

    // Set the CCDNR keyword.
    fits_update_key(fptr, TINT, "CCDNR", &par.CCDNr, "", &status);
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

      // Separate the time value into ots (on-board second clock)
      // and fracsec (subsecond clock) integer values.
      long ots=(long)event.time;
      long fracsec=(long)((event.time-ots)*1.e6); // [micro seconds]

      fits_write_col(fptr, TLONG, cots, row+1, 1, 1, &ots, &status);
      fits_write_col(fptr, TLONG, cfracsec, row+1, 1, 1, &fracsec, &status);
      fits_write_col(fptr, TLONG, cframe, row+1, 1, 1, &event.frame, &status);
      fits_write_col(fptr, TLONG, cpha, row+1, 1, 1, &event.pha, &status);
      int rawx=event.rawx+1;
      fits_write_col(fptr, TINT, crawx, row+1, 1, 1, &rawx, &status);
      int rawy=event.rawy+1;
      fits_write_col(fptr, TINT, crawy, row+1, 1, 1, &rawy, &status);

      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);
    // END of loop over all events in the input file.

    // Append a check sum to the FITS header of the event extension.
    int hdutype;
    fits_movabs_hdu(fptr, 2, &hdutype, &status);
    fits_write_chksum(fptr, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Close the files.
  freeEventFile(&elf, &status);
  if (NULL!=fptr) fits_close_file(fptr, &status);

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


int getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // check if any obsolete keywords are given
  sixt_check_obsolete_keyword(&status);
  CHECK_STATUS_RET(status,EXIT_FAILURE);


  status=ape_trad_query_file_name("RawData", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the input file");
    return(status);
  }
  strcpy(par->RawData, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("eroEvtFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the output file");
    return(status);
  }
  strcpy(par->eroEvtFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_int("CCDNr", &par->CCDNr);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the CCDNr parameter");
    return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  return(status);
}
