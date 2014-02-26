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

#include "makelc.h"


////////////////////////////////////
/** Main procedure. */
int makelc_main() {
  // Program parameters.
  struct Parameters par; 

  // Input event file.
  fitsfile* infptr=NULL;

  // Output light curve.
  long* counts=NULL;
  fitsfile* outfptr=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL:
  set_toolname("makelc");
  set_toolversion("0.09");


  do {  // Beginning of the ERROR handling loop.

    // --- Initialization ---

    // Read the program parameters using PIL library.
    status=makelc_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Set the input event file.
    fits_open_table(&infptr, par.EventList, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Determine timing keywords.
    char comment[MAXMSG];
    char telescop[MAXMSG];
    char instrume[MAXMSG];
    char filter[MAXMSG];
    double mjdref, timezero;
    fits_read_key(infptr, TSTRING, "TELESCOP", telescop, comment, &status);
    fits_read_key(infptr, TSTRING, "INSTRUME", instrume, comment, &status);
    fits_read_key(infptr, TSTRING, "FILTER", filter, comment, &status);
    fits_read_key(infptr, TDOUBLE, "MJDREF", &mjdref, comment, &status);
    fits_read_key(infptr, TDOUBLE, "TIMEZERO", &timezero, comment, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the column containing the time information.
    int ctime;
    fits_get_colnum(infptr, CASEINSEN, "TIME", &ctime, &status);
    CHECK_STATUS_BREAK(status);

    // If the only events within a certain energy band should be
    // considered, determine the column containing the energy /
    // signal information.
    int csignal=0;
    if ((par.Emin>=0.)||(par.Emax>=0.)) {
      fits_write_errmark();
      int opt_status=EXIT_SUCCESS;
      fits_get_colnum(infptr, CASEINSEN, "ENERGY", &csignal, &opt_status);
      fits_clear_errmark();
      if (EXIT_SUCCESS!=opt_status) {
	fits_get_colnum(infptr, CASEINSEN, "SIGNAL", &csignal, &status);
	if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("could not find column 'ENERGY'/'SIGNAL'");
	  break;
	}
      }
    }

    // If the only events within a certain channel range should be
    // considered, determine the column containing the energy /
    // signal information.
    int cpha=0;
    if ((par.Chanmin>=0)||(par.Chanmax>=0)) {
      fits_write_errmark();
      int opt_status=EXIT_SUCCESS;
      fits_get_colnum(infptr, CASEINSEN, "PI", &cpha, &opt_status);
      fits_clear_errmark();
      if (EXIT_SUCCESS!=opt_status) {
	fits_get_colnum(infptr, CASEINSEN, "PHA", &cpha, &status);
	if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("could not find column 'PI'/'PHA'");
	  break;
	}
      }
    }

    // Determine the number of rows in the input file.
    long nrows;
    fits_get_num_rows(infptr, &nrows, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the number of bins for the light curve.
    long nbins=(long)(par.length/par.dt);

    // Allocate memory for the count histogram.
    headas_chat(5, "create light curve with %ld bins ...\n",
		nbins);
    counts=(long*)malloc(nbins*sizeof(long));
    CHECK_NULL_BREAK(counts, status, 
		     "memory allocation for light curve failed");

    // Initialize the light curve with 0.
    long ii; 
    for (ii=0; ii<nbins; ii++) {
      counts[ii]=0;
    }

    // --- END of Initialization ---


    // --- Begin Spectrum Binning ---
    headas_chat(3, "calculate light curve ...\n");

    // LOOP over all events in the FITS table.
    for (ii=0; ii<nrows; ii++) {
      
      // Read the time of the next event from the file.
      double time;
      double dnull=0.0;
      int anynul=0;
      fits_read_col(infptr, TDOUBLE, ctime, ii+1, 1, 1, 
		    &dnull, &time, &anynul, &status);
      CHECK_STATUS_BREAK(status);

      // If the event was detected before the start of the light
      // curve, we have to neglect it.
      if (time<par.TSTART) continue;
      if (time>par.TSTART+par.length) continue;

      // If necessary, read the energy/signal of the event.
      if (csignal>0) {
	float signal;
	float fnull=0.0;
	fits_read_col(infptr, TFLOAT, csignal, ii+1, 1, 1, 
		      &fnull, &signal, &anynul, &status);
	CHECK_STATUS_BREAK(status);

	// Check if the energy of the event lies within the 
	// requested range.
	if ((signal<par.Emin)||(signal>par.Emax)) continue;
      }

      // If necessary, read the energy/signal of the event.
      if (cpha>0) {
	long pha;
	long lnull=0;
	fits_read_col(infptr, TLONG, cpha, ii+1, 1, 1, 
		      &lnull, &pha, &anynul, &status);
	CHECK_STATUS_BREAK(status);

	// Check if the energy of the event lies within the 
	// requested range.
	if ((pha<par.Chanmin)||(pha>par.Chanmax)) continue;
      }
      
      // Determine the respective bin in the light curve.
      long bin=((long)((time-par.TSTART)/par.dt+1.0))-1;
		
      // If the event exceeds the end of the light curve, neglect it.
      if (bin>=nbins) continue;
      
      // Add the event to the light curve.
      assert(bin>=0);
      counts[bin]++;
    }
    CHECK_STATUS_BREAK(status);
    // END of loop over all events in the input file.

    // Store the light curve in the output file.
    headas_chat(3, "store light curve ...\n");

    // Check if the file already exists.
    int exists;
    fits_file_exists(par.LightCurve, &exists, &status);
    CHECK_STATUS_BREAK(status);
    if (0!=exists) {
      if (0!=par.clobber) {
	// Delete the file.
	remove(par.LightCurve);
      } else {
	// Throw an error.
	char msg[MAXMSG];
	sprintf(msg, "file '%s' already exists", par.LightCurve);
	SIXT_ERROR(msg);
	status=EXIT_FAILURE;
	break;
      }
    }

    // Create a new FITS-file.
    char buffer[MAXFILENAME];
    sprintf(buffer, "%s(%s%s)", par.LightCurve, SIXT_DATA_PATH, 
	    "/templates/makelc.tpl");
    fits_create_file(&outfptr, buffer, &status);
    CHECK_STATUS_BREAK(status);

    // Move to the HDU containing the binary table.
    int hdutype;
    fits_movabs_hdu(outfptr, 2, &hdutype, &status);
    CHECK_STATUS_BREAK(status);

    // Get column numbers.
    int ccounts;
    fits_get_colnum(outfptr, CASEINSEN, "COUNTS", &ccounts, &status);
    CHECK_STATUS_BREAK(status);

    // Write header keywords.
    fits_update_key(outfptr, TSTRING, "TELESCOP", telescop,
		    "Telescope name", &status);
    fits_update_key(outfptr, TSTRING, "INSTRUME", instrume,
		    "Instrument name", &status);
    fits_update_key(outfptr, TSTRING, "FILTER", filter,
		    "Filter used", &status);
    fits_update_key(outfptr, TSTRING, "TIMEUNIT", "s",
		    "time unit", &status);
    fits_update_key(outfptr, TDOUBLE, "TIMEDEL", &par.dt,
		    "time resolution", &status);
    fits_update_key(outfptr, TDOUBLE, "MJDREF", &mjdref,
		    "reference MJD", &status);
    fits_update_key(outfptr, TDOUBLE, "TIMEZERO", &timezero,
		    "time offset", &status);
    float timepixr=0.f;
    fits_update_key(outfptr, TFLOAT, "TIMEPIXR", &timepixr,
		    "time stamp at beginning of bin", &status);
    fits_update_key(outfptr, TDOUBLE, "TSTART", &par.TSTART,
		    "start time", &status);
    double dbuffer=par.TSTART+par.length;
    fits_update_key(outfptr, TDOUBLE, "TSTOP", &dbuffer,
		    "stop time", &status);
    fits_update_key(outfptr, TFLOAT, "E_MIN", &par.Emin,
		    "low energy for channel (keV)", &status);
    fits_update_key(outfptr, TFLOAT, "E_MAX", &par.Emax,
		    "high energy for channel (keV)", &status);
    CHECK_STATUS_BREAK(status);

    // The ouput table does not contain a TIME column. The 
    // beginning (TIMEPIXR=0.0) of the n-th time bin (n>=1) 
    // is determined as t(n)=TIMEZERO + TIMEDEL*(n-1).

    // Write the data into the table.
    for (ii=0; ii<nbins; ii++) {
      fits_write_col(outfptr, TLONG, ccounts, ii+1, 1, 1, 
		     &(counts[ii]), &status);
      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Close the files.
  if (NULL!= infptr) fits_close_file( infptr, &status);
  if (NULL!=outfptr) fits_close_file(outfptr, &status);

  // Release memory.
  if (NULL!=counts) free(counts);

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


int makelc_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("EventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the event list file");
    return(status);
  } 
  strcpy(par->EventList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("LightCurve", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the output light curve file");
    return(status);
  } 
  strcpy(par->LightCurve, sbuffer);
  free(sbuffer);

  status=ape_trad_query_double("TSTART", &par->TSTART);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the start time of the light curve");
    return(status);
  } 

  status=ape_trad_query_double("length", &par->length);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the length of the light curve");
    return(status);
  } 

  status=ape_trad_query_double("dt", &par->dt);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the time resolution of the light curve");
    return(status);
  } 

  status=ape_trad_query_float("Emin", &par->Emin);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the lower boundary of the energy band");
    return(status);
  } 

  status=ape_trad_query_float("Emax", &par->Emax);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the upper boundary of the energy band");
    return(status);
  } 

  status=ape_trad_query_long("Chanmin", &par->Chanmin);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the lower boundary of the channel range");
    return(status);
  } 

  status=ape_trad_query_long("Chanmax", &par->Chanmax);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the upper boundary of the channel range");
    return(status);
  } 

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  return(status);
}

