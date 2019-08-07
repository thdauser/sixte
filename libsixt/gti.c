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

#include "gti.h"


GTI* newGTI(int* const status)
{
  GTI* gti=(GTI*)malloc(sizeof(GTI));
  CHECK_NULL_RET(gti, *status,
		 "memory allocation for GTI data structure failed", gti);

  // Initialize.
  HDgti_init(gti);

  return(gti);
}


void freeGTI(GTI** const gti)
{
  if (NULL!=*gti) {
    HDgti_free(*gti);
    free(*gti);
    *gti=NULL;
  }
}


GTI* loadGTI(char* const filename, int* const status)
{
  GTI* gti=NULL;

  // Allocate memory.
  gti=newGTI(status);
  CHECK_NULL(gti, *status,
	     "memory allocation for GTI data structure failed");

  // Try to load an extension with the name 'GTI' or 'STDGTI'.
  int status2=EXIT_SUCCESS;
  fits_write_errmark();
  HDgti_read(filename, gti, "GTI", "START", "STOP", NULL, NULL, &status2);
  if (EXIT_SUCCESS!=status2) {
    status2=EXIT_SUCCESS;
    HDgti_read(filename, gti, "STDGTI", "START", "STOP", NULL, NULL, &status2);
  }
  fits_clear_errmark();
  if (EXIT_SUCCESS!=status2) {
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "could not find extension 'GTI' or 'STDGTI' "
	    "in file '%s'", filename);
    SIXT_ERROR(msg);
    return(NULL);
  }

  // Make sure that TIMEZERO==0.0.
  verifyTIMEZERO(gti->timezero, status);
  CHECK_STATUS_RET(*status, gti);

  return(gti);
}


void saveGTI(GTI* const gti,
	     const char* const filename,
	     const char clobber,
	     int* const status)
{
  fitsfile* fptr=NULL;

  // Check if the file already exists.
  int exists;
  fits_file_exists(filename, &exists, status);
  CHECK_STATUS_VOID(*status);
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
      return;
    }
  }

  // Create a new FITS file.
  fits_create_file(&fptr, filename, status);
  CHECK_STATUS_VOID(*status);

  // Store the GTI extension.
  saveGTIExt(fptr, "GTI", gti, status);
  CHECK_STATUS_VOID(*status);

  // Set the time-keyword in the event list header.
  char datestr[MAXMSG];
  int timeref;
  fits_get_system_time(datestr, &timeref, status);
  CHECK_STATUS_VOID(*status);
  fits_update_key(fptr, TSTRING, "DATE", datestr,
		  "File creation date", status);
  CHECK_STATUS_VOID(*status);

  // Close the FITS file.
  fits_close_file(fptr, status);
  CHECK_STATUS_VOID(*status);
}


void saveGTIExt(fitsfile* const fptr,
		char* const extname,
		GTI* const gti,
		int* const status)
{
	char buffer_telescop[MAXMSG];
	char buffer_instrume[MAXMSG];


	// Check if the EVENTS extension already exists and move fptr to that extension
	fits_movnam_hdu(fptr, BINARY_TBL, "EVENTS", 0, status);
	if (*status != EXIT_SUCCESS) {
		char msg[MAXMSG];
		sprintf(msg, "Unable to write keywords INSTRUME, TELESCOP to %s extension!",extname);
		SIXT_WARNING(msg);

		fits_clear_errmsg();
		*status=EXIT_SUCCESS;
	}
	else{
		// Read INSTRUME and TELESCOP keywords from EVENTS extension
		fits_read_key(fptr, TSTRING, "INSTRUME", buffer_instrume, NULL, status);
		fits_read_key(fptr, TSTRING, "TELESCOP", buffer_telescop, NULL, status);
	}

	HDgti_write(fptr, gti, extname, "START", "STOP", status);

	if( strlen(buffer_telescop) && strlen(buffer_instrume) ){
		// Move to STDGTI extension
		fits_movnam_hdu(fptr, BINARY_TBL, extname, 0, status);
		if (*status != EXIT_SUCCESS) {
			char buffer[MAXMSG];
			fits_get_errstatus(*status, buffer);
			SIXT_WARNING(buffer);
			CHECK_STATUS_VOID(*status);
		}

		// Write TELESCOP and INSTRUME keywords to STDGTI extension
		fits_update_key(fptr, TSTRING, "INSTRUME", buffer_instrume, NULL, status);
		fits_update_key(fptr, TSTRING, "TELESCOP", buffer_telescop, NULL, status);
	}
}


void appendGTI(GTI* const gti,
	       const double start,
	       const double stop,
	       int* const status)
{
  // Allocate memory.
  if (gti->ngti>=gti->maxgti) {
    HDgti_grow(gti, gti->ngti+1, status);
    CHECK_STATUS_VOID(*status);
  }

  // Store the start and the stop time.
  gti->start[gti->ngti]=start;
  gti->stop[gti->ngti] =stop;
  gti->ngti++;
}


double sumGTI(GTI* const gti)
{
  double sum=0.;
  int ii;
  for (ii=0; ii<gti->ngti; ii++) {
    sum+=gti->stop[ii]-gti->start[ii];
  }
  return(sum);
}


GTI* getGTIFromFileOrContinuous(char* const filename,
				const double tstart,
				const double tstop,
				const double mjdref,
				int* const status)
{
  GTI* gti=NULL;

  // If available, load the specified GTI file.
  if (strlen(filename)>0) {
    char ucase_buffer[MAXFILENAME];
    strcpy(ucase_buffer, filename);
    strtoupper(ucase_buffer);
    if (0!=strcmp(ucase_buffer, "NONE")) {
      gti=loadGTI(filename, status);
      CHECK_STATUS_RET(*status, gti);
      verifyMJDREF(mjdref, gti->mjdref, "in GTI file", status);
      CHECK_STATUS_RET(*status, gti);
      SIXT_WARNING("the given GTI file overwrites the settings for TSTART, MJDREF and Exposure");
      printf("    -> simulating from TSTART=%.6e until TSTOP=%.6e (MJDREF=%.3f) \n",
    		  gti->start[0],gti->stop[gti->ngti-1],gti->mjdref);
    }
  }

  // If not, create a dummy GTI from TSTART and TSTOP.
  if (NULL==gti) {
    gti=newGTI(status);
    CHECK_STATUS_RET(*status, gti);
    gti->mjdref=mjdref;
    gti->timezero=0.0;
    appendGTI(gti, tstart, tstop, status);
    CHECK_STATUS_RET(*status, gti);
  }

  return(gti);
}
