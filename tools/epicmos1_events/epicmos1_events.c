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


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

#include "epicmos1_events.h"


int epicmos1_events_main()
{
  // Containing all programm parameters read by PIL.
  struct Parameters par;

  // Input event file.
  EventFile* elf=NULL;

  // File pointer to the output EPIC-mos1 event file.
  fitsfile* fptr=NULL;

  // WCS data structure used for projection.
  struct wcsprm wcs={ .flag=-1 };
  // String buffer for FITS header.
  char* headerstr=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  // Plate scale [0.05arcsec/pixel].
  const double ps=atan(40.e-6/7.5)*180./M_PI*3600./0.05;

  // Register HEATOOL:
  set_toolname("epicmos1_events");
  set_toolversion("0.01");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    headas_chat(3, "initialization ...\n");

    // Read parameters using PIL library:
    if ((status=getpar(&par))) break;

    // Open the input event file.
    elf=openEventFile(par.EvtFile, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Check if the input file contains recombined event patterns.
    char evtype[MAXMSG], comment[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "EVTYPE", evtype, comment, &status);
    if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("could not read FITS keyword 'EVTYPE'");
      break;
    }
    strtoupper(evtype);
    if (0!=strcmp(evtype, "PATTERN")) {
      status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "event type of input file is '%s' (must be 'PATTERN')", evtype);
      SIXT_ERROR(msg);
      break;
    }

    // Read keywords from the input file.
    double mjdref=0.0;
    fits_read_key(elf->fptr, TDOUBLE, "MJDREF", &mjdref, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'MJDREF' from input "
	      "event list '%s'", par.EvtFile);
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
	      "event list '%s'", par.EvtFile);
      SIXT_ERROR(msg);
      break;
    }

    char time_obs[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "TIME-OBS", time_obs, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'TIME-OBS' from input "
	      "event list '%s'", par.EvtFile);
      SIXT_ERROR(msg);
      break;
    }

    char date_end[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "DATE-END", date_end, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'DATE-END' from input "
	      "event list '%s'", par.EvtFile);
      SIXT_ERROR(msg);
      break;
    }

    char time_end[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "TIME-END", time_end, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'TIME-END' from input "
	      "event list '%s'", par.EvtFile);
      SIXT_ERROR(msg);
      break;
    }

    double tstart=0.0;
    fits_read_key(elf->fptr, TDOUBLE, "TSTART", &tstart, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'TSTART' from input "
	      "event list '%s'", par.EvtFile);
      SIXT_ERROR(msg);
      break;
    }

    double tstop=0.0;
    fits_read_key(elf->fptr, TDOUBLE, "TSTOP", &tstop, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'TSTOP' from input "
	      "event list '%s'", par.EvtFile);
      SIXT_ERROR(msg);
      break;
    }

    char ancrfile[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "ANCRFILE", ancrfile, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'ANCRFILE' from input "
	      "event list '%s'", par.EvtFile);
      SIXT_ERROR(msg);
      break;
    }

    char respfile[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "RESPFILE", respfile, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'RESPFILE' from input "
	      "event list '%s'", par.EvtFile);
      SIXT_ERROR(msg);
      break;
    }

    float ra_pnt=0.0;
    fits_read_key(elf->fptr, TFLOAT, "RA_PNT", &ra_pnt, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'RA_PNT' from input "
	      "event list '%s'", par.EvtFile);
      SIXT_ERROR(msg);
      break;
    }

    float dec_pnt=0.0;
    fits_read_key(elf->fptr, TFLOAT, "DEC_PNT", &dec_pnt, comment, &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not read FITS keyword 'DEC_PNT' from input "
	      "event list '%s'", par.EvtFile);
      SIXT_ERROR(msg);
      break;
    }

    // Verify values of MJDREF and TIMEZERO.
    verifyMJDREF(xmmmjdref, mjdref, "in event file", &status);
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
    fits_file_exists(par.EPICmos1EventList, &exists, &status);
    CHECK_STATUS_BREAK(status);
    if (0!=exists) {
      if (0!=par.clobber) {
	// Delete the file.
	remove(par.EPICmos1EventList);
      } else {
	// Throw an error.
	char msg[MAXMSG];
	sprintf(msg, "file '%s' already exists", par.EPICmos1EventList);
	SIXT_ERROR(msg);
	status=EXIT_FAILURE;
	break;
      }
    }

    // Create and open a new FITS file.
    headas_chat(3, "create new EPIC-mos1 event list file '%s' ...\n",
		par.EPICmos1EventList);
    fits_create_file(&fptr, par.EPICmos1EventList, &status);
    CHECK_STATUS_BREAK(status);

    // Create the event table.
    char *ttype[]={"TIME",
		   "RAWX", "RAWY",
		   "DETX", "DETY", "X", "Y",
		   "PHA", "PI", "FLAG", "PATTERN"};
    char *tform[]={"D",
		   "I", "I",
		   "I", "I", "J", "J",
		   "I", "I", "J", "B"};
    char *tunit[]={"s",
		   "pixel", "pixel",
		   "0.05arcsec", "0.05arcsec", "0.05arcsec", "0.05arcsec",
		   "channel", "eV", "", ""};
    fits_create_tbl(fptr, BINARY_TBL, 0, 11, ttype, tform, tunit,
		    "EVENTS", &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not create binary table for events "
	      "in file '%s'", par.EPICmos1EventList);
      SIXT_ERROR(msg);
      break;
    }

    // Insert the standard header keywords.
    sixt_add_fits_stdkeywords_obsolete(fptr, 1, "XMM", "EM1", "", ancrfile, respfile,
			      mjdref, timezero, tstart, tstop, &status);
    CHECK_STATUS_BREAK(status);
    sixt_add_fits_stdkeywords_obsolete(fptr, 2, "XMM", "EM1", "", ancrfile, respfile,
			      mjdref, timezero, tstart, tstop, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the column numbers.
    int ctime, crawx, crawy, cdetx, cdety, cx, cy,
      cpha, cpi, cflag, cpattern;
    fits_get_colnum(fptr, CASEINSEN, "TIME", &ctime, &status);
    fits_get_colnum(fptr, CASEINSEN, "RAWX", &crawx, &status);
    fits_get_colnum(fptr, CASEINSEN, "RAWY", &crawy, &status);
    fits_get_colnum(fptr, CASEINSEN, "DETX", &cdetx, &status);
    fits_get_colnum(fptr, CASEINSEN, "DETY", &cdety, &status);
    fits_get_colnum(fptr, CASEINSEN, "X", &cx, &status);
    fits_get_colnum(fptr, CASEINSEN, "Y", &cy, &status);
    fits_get_colnum(fptr, CASEINSEN, "PHA", &cpha, &status);
    fits_get_colnum(fptr, CASEINSEN, "PI", &cpi, &status);
    fits_get_colnum(fptr, CASEINSEN, "FLAG", &cflag, &status);
    fits_get_colnum(fptr, CASEINSEN, "PATTERN", &cpattern, &status);
    CHECK_STATUS_BREAK(status);

    // Set the TLMIN and TLMAX keywords.
    char keyword[MAXMSG];
    // For the columns RAWX and RAWY.
    long tlmin_rawx=1, tlmax_rawx=384;
    long tlmin_rawy=1, tlmax_rawy=400;
    sprintf(keyword, "TLMIN%d", crawx);
    fits_update_key(fptr, TLONG, keyword, &tlmin_rawx, "", &status);
    sprintf(keyword, "TLMAX%d", crawx);
    fits_update_key(fptr, TLONG, keyword, &tlmax_rawx, "", &status);
    sprintf(keyword, "TLMIN%d", crawy);
    fits_update_key(fptr, TLONG, keyword, &tlmin_rawy, "", &status);
    sprintf(keyword, "TLMAX%d", crawy);
    fits_update_key(fptr, TLONG, keyword, &tlmax_rawy, "", &status);
    CHECK_STATUS_BREAK(status);

    // TODO DETX, DETY, X, Y.

    // For the column PHA.
    int tlmin_pha=0, tlmax_pha=4095;
    sprintf(keyword, "TLMIN%d", cpha);
    fits_update_key(fptr, TINT, keyword, &tlmin_pha, "", &status);
    sprintf(keyword, "TLMAX%d", cpha);
    fits_update_key(fptr, TINT, keyword, &tlmax_pha, "", &status);
    CHECK_STATUS_BREAK(status);

    // For the column PI.
    int tlmin_pi=0, tlmax_pi=20480;
    sprintf(keyword, "TLMIN%d", cpi);
    fits_update_key(fptr, TINT, keyword, &tlmin_pi, "", &status);
    sprintf(keyword, "TLMAX%d", cpi);
    fits_update_key(fptr, TINT, keyword, &tlmax_pi, "", &status);
    CHECK_STATUS_BREAK(status);


    // Set up the WCS data structure.
    if (0!=wcsini(1, 2, &wcs)) {
      SIXT_ERROR("initalization of WCS data structure failed");
      status=EXIT_FAILURE;
      break;
    }
    wcs.naxis=2;
    wcs.crpix[0]=0.0;
    wcs.crpix[1]=0.0;
    wcs.crval[0]=ra_pnt;
    wcs.crval[1]=dec_pnt;
    wcs.cdelt[0]=-0.05/3600.;
    wcs.cdelt[1]= 0.05/3600.;
    strcpy(wcs.cunit[0], "deg");
    strcpy(wcs.cunit[1], "deg");
    strcpy(wcs.ctype[0], "RA---TAN");
    strcpy(wcs.ctype[1], "DEC--TAN");

    // Update the WCS keywords in the output file.
    sprintf(keyword, "TCTYP%d", cx);
    fits_update_key(fptr, TSTRING, keyword, wcs.ctype[0],
		    "projection type", &status);
    sprintf(keyword, "TCTYP%d", cy);
    fits_update_key(fptr, TSTRING, keyword, wcs.ctype[1],
		    "projection type", &status);
    sprintf(keyword, "TCRVL%d", cx);
    fits_update_key(fptr, TDOUBLE, keyword, &wcs.crval[0],
		    "reference value", &status);
    sprintf(keyword, "TCRVL%d", cy);
    fits_update_key(fptr, TDOUBLE, keyword, &wcs.crval[1],
		    "reference value", &status);
    sprintf(keyword, "TCRPX%d", cx);
    fits_update_key(fptr, TFLOAT, keyword, &wcs.crpix[0],
		    "reference point", &status);
    sprintf(keyword, "TCRPX%d", cy);
    fits_update_key(fptr, TFLOAT, keyword, &wcs.crpix[1],
		    "reference point", &status);
    sprintf(keyword, "TCDLT%d", cx);
    fits_update_key(fptr, TDOUBLE, keyword, &wcs.cdelt[0],
		    "pixel increment", &status);
    sprintf(keyword, "TCDLT%d", cy);
    fits_update_key(fptr, TDOUBLE, keyword, &wcs.cdelt[1],
		    "pixel increment", &status);
    sprintf(keyword, "TCUNI%d", cx);
    fits_update_key(fptr, TSTRING, keyword, wcs.cunit[0],
		    "axis units", &status);
    sprintf(keyword, "TCUNI%d", cy);
    fits_update_key(fptr, TSTRING, keyword, wcs.cunit[1],
		    "axis units", &status);
    CHECK_STATUS_BREAK(status);

    fits_update_key(fptr, TSTRING, "REFXCTYP", wcs.ctype[0],
		    "projection type", &status);
    fits_update_key(fptr, TSTRING, "REFYCTYP", wcs.ctype[1],
		    "projection type", &status);
    fits_update_key(fptr, TSTRING, "REFXCUNI", wcs.cunit[0],
		    "axis units", &status);
    fits_update_key(fptr, TSTRING, "REFYCUNI", wcs.cunit[1],
		    "axis units", &status);
    fits_update_key(fptr, TFLOAT, "REFXCRPX", &wcs.crpix[0],
		    "reference value", &status);
    fits_update_key(fptr, TFLOAT, "REFYCRPX", &wcs.crpix[1],
		    "reference value", &status);
    fits_update_key(fptr, TDOUBLE, "REFXCRVL", &wcs.crval[0],
		    "reference value", &status);
    fits_update_key(fptr, TDOUBLE, "REFYCRVL", &wcs.crval[1],
		    "reference value", &status);
    fits_update_key(fptr, TDOUBLE, "REFXCDLT", &wcs.cdelt[0],
		    "pixel increment", &status);
    fits_update_key(fptr, TDOUBLE, "REFYCDLT", &wcs.cdelt[1],
		    "pixel increment", &status);
    CHECK_STATUS_BREAK(status);

    // --- END of initialization ---

    // --- Beginning of copy events ---

    headas_chat(3, "copy events ...\n");

    // Actual minimum and maximum values of X and Y.
    long refxdmin, refxdmax, refydmin, refydmax;

    // Loop over all events in the FITS file.
    long row;
    for (row=0; row<elf->nrows; row++) {

      // Read the next event from the input file.
      Event event;
      getEventFromFile(elf, row+1, &event, &status);
      CHECK_STATUS_BREAK(status);

      // Determine the event data based on the event information.
      EPICmos1Event ev;

      // Time.
      ev.time=event.time;

      // Convert world coordinates to image coordinates X and Y.
      double world[2]={ event.ra*180./M_PI, event.dec*180./M_PI };
      double imgcrd[2], pixcrd[2];
      double phi, theta;
      int wcsstatus=0;
      wcss2p(&wcs, 1, 2, world, &phi, &theta, imgcrd, pixcrd, &wcsstatus);
      if (0!=wcsstatus) {
	char msg[MAXMSG];
	sprintf(msg,
		"WCS coordinate conversion failed (RA=%lf, Dec=%lf, error code %d)",
		world[0], world[1], wcsstatus);
	SIXT_ERROR(msg);
	status=EXIT_FAILURE;
	break;
      }
      ev.x=(long)pixcrd[0];
      if (pixcrd[0] < 0.) ev.x--;
      ev.y=(long)pixcrd[1];
      if (pixcrd[1] < 0.) ev.y--;

      // Determine the actual minimum and maximum values of X and Y.
      if (0==row) {
	refxdmin=ev.x;
	refxdmax=ev.x;
	refydmin=ev.y;
	refydmax=ev.y;
      }
      if (ev.x<refxdmin) {
	refxdmin=ev.x;
      }
      if (ev.x>refxdmax) {
	refxdmax=ev.x;
      }
      if (ev.y<refydmin) {
	refydmin=ev.y;
      }
      if (ev.y>refydmax) {
	refydmax=ev.y;
      }

      // Pixel coordinates.
      ev.rawx=event.rawx+1;
      ev.rawy=event.rawy+1;
      ev.detx=(event.rawx+1-160.5)*ps;
      ev.dety=(event.rawy+1-220.5)*ps;

      // Pattern.
      ev.pattern=event.type;

      // TODO In the current implementation the value of FLAG is only
      // set to 0 or 2. This needs to be extended later.
      if (event.type>=0) {
	ev.flag=0;
      } else {
	// Invalid pattern.
	ev.flag=2;
      }

      // PI and PHA.
      ev.pha=event.pha;
      ev.pi =(int)(event.signal*1000.);

      // Store the event in the output file.
      fits_write_col(fptr, TDOUBLE, ctime, row+1, 1, 1, &ev.time, &status);
      fits_write_col(fptr, TINT, crawx, row+1, 1, 1, &ev.rawx, &status);
      fits_write_col(fptr, TINT, crawy, row+1, 1, 1, &ev.rawy, &status);
      fits_write_col(fptr, TINT, cdetx, row+1, 1, 1, &ev.detx, &status);
      fits_write_col(fptr, TINT, cdety, row+1, 1, 1, &ev.dety, &status);
      fits_write_col(fptr, TLONG, cx, row+1, 1, 1, &ev.x, &status);
      fits_write_col(fptr, TLONG, cy, row+1, 1, 1, &ev.y, &status);
      fits_write_col(fptr, TINT, cpha, row+1, 1, 1, &ev.pha, &status);
      fits_write_col(fptr, TINT, cpi, row+1, 1, 1, &ev.pi, &status);
      fits_write_col(fptr, TLONG, cflag, row+1, 1, 1, &ev.flag, &status);
      fits_write_col(fptr, TBYTE, cpattern, row+1, 1, 1, &ev.pattern, &status);
      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);
    // END of loop over all events in the FITS file.

    // Set the REF?DMIN/MAX keywords.
    fits_update_key(fptr, TLONG, "REFXDMIN", &refxdmin, "", &status);
    fits_update_key(fptr, TLONG, "REFXDMAX", &refxdmax, "", &status);
    fits_update_key(fptr, TLONG, "REFYDMIN", &refydmin, "", &status);
    fits_update_key(fptr, TLONG, "REFYDMAX", &refydmax, "", &status);
    CHECK_STATUS_BREAK(status);

    // --- End of copy events ---

    // Append a check sum to the header of the event extension.
    int hdutype=0;
    fits_movabs_hdu(fptr, 2, &hdutype, &status);
    fits_write_chksum(fptr, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Close the files.
  freeEventFile(&elf, &status);
  if (NULL!=fptr) fits_close_file(fptr, &status);

  // Release memory.
  wcsfree(&wcs);
  if (NULL!=headerstr) free(headerstr);

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

  status=ape_trad_query_file_name("EvtFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the input pattern list");
    return(status);
  }
  strcpy(par->EvtFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("EPICmos1EventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the output event list");
    return(status);
  }
  strcpy(par->EPICmos1EventList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  return(status);
}
