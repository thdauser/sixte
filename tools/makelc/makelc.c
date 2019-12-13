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

#include "makelc.h"



////////////////////////////////////
/** Main procedure. */
int makelc_main() {
  // Program parameters.
  struct Parameters par;

  // Input event file.
  fitsfile* infptr=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  // Register HEATOOL:
  set_toolname("makelc");
  set_toolversion("0.10");

  do {  // Beginning of the ERROR handling loop.

    // --- Initialization ---

    // Read the program parameters using PIL library.
    status=makelc_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Set the input event file.
    fits_open_table(&infptr, par.EvtFile, READONLY, &status);
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

    
    // If a GTI file is given, overwrite TSTART and length, and write
    // multiple output lightcurves.
    long nrows_gti; // No. rows in GTI file, if none given set to 1 later
    double* starts=NULL; // START values from GTI, if none given set to par.TSTART
    double* lengths=NULL; // STOP-START values from GTI, if non fiven set to par.length

    // Get list of START, STOP values from GTI file and compute length
    if (0!=strcmp(par.GTIFile,"none")) {
      fitsfile *gtifptr=NULL;

      // Open the GTI file.
      fits_open_file(&gtifptr, par.GTIFile, READONLY, &status);
      fits_movnam_hdu(gtifptr, BINARY_TBL, "GTI", 0, &status);
      CHECK_STATUS_BREAK(status);

      // Determine the column containing the START and STOP (s) information.
      int cstart; int cstop;
      fits_get_colnum(gtifptr, CASEINSEN, "START", &cstart, &status);
      fits_get_colnum(gtifptr, CASEINSEN, "STOP", &cstop, &status);
      CHECK_STATUS_BREAK(status);

      // Determine the number of rows in the GTI file.
      fits_get_num_rows(gtifptr, &nrows_gti, &status);
      CHECK_STATUS_BREAK(status);
      
      // Initialize memory for the time arrays of the GTI file.
      starts=(double*)malloc(nrows_gti*sizeof(double));
      lengths=(double*)malloc(nrows_gti*sizeof(double));
      CHECK_NULL_BREAK(starts, status,
		       "memory allocation for start array in GTI reading failed");
      CHECK_NULL_BREAK(lengths, status,
		       "memory allocation for lengths array in GTI reading failed");
      
      // LOOP over all start and stop times in the GTI FITS table.
      int ii;
      for (ii=0; ii<nrows_gti; ii++) {
	double dnull=0.0;
	int anynul=0;
	double stop;
	fits_read_col(gtifptr, TDOUBLE, cstart, ii+1, 1, 1,
		      &dnull, &starts[ii], &anynul, &status);
	fits_read_col(gtifptr, TDOUBLE, cstop, ii+1, 1, 1,
		      &dnull, &stop, &anynul, &status);
	CHECK_STATUS_BREAK(status);
	lengths[ii]=stop-starts[ii];
      }
      // Close GTI file
      if (NULL!=gtifptr) fits_close_file(gtifptr, &status);
    }
    
    
    // If no GTI file given take TSTART and length --> for loop
    // executed only once.
    if (0==strcmp(par.GTIFile,"none")) {
      nrows_gti=1;
      // Initialize memory for the time arrays of length 1
      starts=(double*)malloc(sizeof(double));
      lengths=(double*)malloc(sizeof(double));
      starts[0]=par.TSTART;
      lengths[0]=par.length;
    } else {
      headas_chat(3, "### Warning: the given GTI file overwrites the settings for TSTART and Length!\n");
    }

    // If a GTI file is provided overwrite TSTART, length and loop
    // over all GTI start and end times and compute lightcurves for
    // each of these intervals.
    int jj;
    // LOOP over all GTI start and length values.
    for (jj=0; jj<nrows_gti; jj++) {

      // Output light curve.
      long* counts=NULL;
      fitsfile* outfptr=NULL;
      
      // Get tstart, length from array (either from GTI or TSTART,length input)
      double tstart=starts[jj];
      double length=lengths[jj];

      // Modify output filename and append _idx in case intervals read from GTI
      char newoutfile[MAXFILENAME];
      if (0==strcmp(par.GTIFile,"none")) {
      	sprintf(newoutfile, "%s", par.LightCurve);
      } else {
	// Get filename without extension
	char *basename; char *extname;
	basename=removeStrExt(par.LightCurve);
	extname=get_filename_ext(par.LightCurve);
	sprintf(newoutfile, "%s_%i.%s", basename, jj, extname);
	free (basename);
      }

      headas_chat(4, "Lightcurve  %i/%i: start=%.2fs, length=%.2fs\n",jj+1,nrows_gti,tstart,length);
      
      // Determine the number of rows in the input file.
      long nrows;
      fits_get_num_rows(infptr, &nrows, &status);
      CHECK_STATUS_BREAK(status);
      
      // Determine the number of bins for the light curve.
      long nbins=(long)(length/par.dt);
      
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
	if (time<tstart) continue;
	if (time>tstart+length) continue;
	
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
	long bin=((long)((time-tstart)/par.dt+1.0))-1;
	
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
      fits_file_exists(newoutfile, &exists, &status);
      CHECK_STATUS_BREAK(status);
      if (0!=exists) {
	if (0!=par.clobber) {
	  // Delete the file.
	  remove(newoutfile);
	} else {
	  // Throw an error.
	  char msg[MAXMSG];
	  sprintf(msg, "file '%s' already exists", newoutfile);
	  SIXT_ERROR(msg);
	  status=EXIT_FAILURE;
	  break;
	}
      }

      // Create a new FITS-file.
      char buffer[MAXFILENAME];
      sprintf(buffer, "%s(%s%s)", newoutfile, SIXT_DATA_PATH,
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
      fits_update_key(outfptr, TDOUBLE, "TSTART", &tstart,
		      "start time", &status);
      double dbuffer=tstart+length;
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

      
      if (NULL!=outfptr) fits_close_file(outfptr, &status);
      // Release memory.
      if (NULL!=counts) free(counts);

      
    } // END of while loop over different GTI times
      
    // --- Cleaning up ---
    headas_chat(3, "cleaning up ...\n");
    
    // Close the files.
    if (NULL!= infptr) fits_close_file( infptr, &status);

    // Release memory.
    if (NULL!=starts) free(starts);
    if (NULL!=lengths) free(lengths);
    
    
  } while(0); // END of the error handling loop.

  
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

  status=ape_trad_query_file_name("EvtFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the event list file");
    return(status);
  }
  strcpy(par->EvtFile, sbuffer);
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

  status=ape_trad_query_file_name("GTIFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the GTI file");
    return(status);
  }
  strcpy(par->GTIFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  return(status);
}


// Function to get filename extension
char *get_filename_ext(char *filename) {
  char *dot = strrchr(filename, '.');
  if(!dot || dot == filename) return "";
  return dot + 1;
}


// Function to get the pathname of string without extension
char *removeStrExt(char* myStr) {
  char *retStr;
  char *lastExt;
  if (myStr == NULL) return NULL;
  if ((retStr = malloc (strlen (myStr) + 1)) == NULL) return NULL;
  strcpy (retStr, myStr);
  lastExt = strrchr (retStr, '.');
  if (lastExt != NULL)
    *lastExt = '\0';
  return retStr;
}
