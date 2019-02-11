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
*/

#include "attitudefile.h"


AttitudeFileEntry read_AttitudeFileEntry(AttitudeFile* const af, int* const status)
{
  AttitudeFileEntry afe={ .time = 0. };
  int anynul=0;
 
  // Time
  fits_read_col(af->fptr, TDOUBLE, af->ctime, af->row+1, 1, 1, 
		&afe.time, &afe.time, &anynul, status);
  // RA
  fits_read_col(af->fptr, TFLOAT, af->cra, af->row+1, 1, 1, 
		&afe.ra, &afe.ra, &anynul, status);
  // Dec
  fits_read_col(af->fptr, TFLOAT, af->cdec, af->row+1, 1, 1, 
		&afe.dec, &afe.dec, &anynul, status);
  // Roll angle.
  if (af->crollang>0) {
    fits_read_col(af->fptr, TFLOAT, af->crollang, af->row+1, 1, 1, 
		  &afe.rollang, &afe.rollang, &anynul, status);
  } else {
    afe.rollang=0.0;
  }

  return(afe);
}


AttitudeFile* open_AttitudeFile(const char filename[], 
				const int access_mode, 
				int* const status)
{
  AttitudeFile* af=NULL;

  do { // Beginning of ERROR handling loop

    af=(AttitudeFile*)malloc(sizeof(AttitudeFile));
    CHECK_NULL_BREAK(af, *status, "memory allocation failed");
    
    // Initialize.
    af->fptr =NULL;
    af->row  =0;
    af->nrows=0;
    af->ctime=0;
    af->cra  =0;
    af->cdec =0;
    af->crollang=0;

    // Open the FITS file table for reading:
    fits_open_table(&af->fptr, filename, access_mode, status);
    CHECK_STATUS_BREAK(*status);

    // Get the HDU type
    int hdutype;
    fits_get_hdu_type(af->fptr, &hdutype, status);
    CHECK_STATUS_BREAK(*status);

    // Image HDU results in an error message.
    if (hdutype==IMAGE_HDU) {
      char msg[MAXMSG];
      *status=EXIT_FAILURE;
      sprintf(msg, "no table extension available in attitude "
	      "file '%s'", filename);
      SIXT_ERROR(msg);
      break;
    }

    // Determine the number of rows in the attitude file.
    fits_get_num_rows(af->fptr, &af->nrows, status);
    CHECK_STATUS_BREAK(*status);


    // Determine the individual column numbers.
    // Required columns.
    fits_get_colnum(af->fptr, CASEINSEN, "TIME", &af->ctime, status);
    CHECK_STATUS_BREAK(*status);

    // Optional columns.
    int opt_status=EXIT_SUCCESS;
    fits_write_errmark();

    fits_get_colnum(af->fptr, CASEINSEN, "RA", &af->cra, &opt_status);
    fits_clear_errmark();
    if (0==af->cra) {
      fits_get_colnum(af->fptr, CASEINSEN, "VIEWRA", &af->cra, status);
      CHECK_STATUS_BREAK(*status);
    }
    opt_status=EXIT_SUCCESS;

    fits_get_colnum(af->fptr, CASEINSEN, "DEC", &af->cdec, &opt_status);
    fits_clear_errmark();
    if (0==af->cdec) {
      fits_get_colnum(af->fptr, CASEINSEN, "VIEWDEC", &af->cdec, status);
      CHECK_STATUS_BREAK(*status);
    }
    opt_status=EXIT_SUCCESS;

    fits_get_colnum(af->fptr, CASEINSEN, "ROLLANG", &af->crollang, &opt_status);
    if (0==af->crollang) {
    	opt_status=EXIT_SUCCESS;
    	fits_get_colnum(af->fptr, CASEINSEN, "ROLL", &af->crollang, &opt_status);
    	if (0==af->crollang) {
    		opt_status=EXIT_SUCCESS;
    		fits_get_colnum(af->fptr, CASEINSEN, "ROLL1", &af->crollang, &opt_status);
    		CHECK_STATUS_BREAK(*status);
    		if (af->crollang==0){
    			headas_chat(3, " *** warning : could not find ROLLANG/ROLL/ROLL1 field in Attitude file %s\n",filename);
    			headas_chat(3, "               setting roll angle to 0.0 %i \n",af->crollang);
    		}
    	}
    }

    fits_clear_errmark();
    // End of determine the column numbers.

  } while(0); // END of ERROR handling loop.

  return(af);
}

