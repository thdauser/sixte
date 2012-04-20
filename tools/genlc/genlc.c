#include "genlc.h"


////////////////////////////////////
/** Main procedure. */
int genlc_main() {
  // Program parameters.
  struct Parameters par; 

  // Input event list file.
  fitsfile* infptr=NULL;

  // Output light curve.
  long* counts=NULL;
  fitsfile* outfptr=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL:
  set_toolname("genlc");
  set_toolversion("0.03");


  do {  // Beginning of the ERROR handling loop.

    // --- Initialization ---

    // Read the program parameters using PIL library.
    status=genlc_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Set the input pattern file.
    fits_open_table(&infptr, par.EventList, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the column containing the time information.
    int citime;
    fits_get_colnum(infptr, CASEINSEN, "TIME", &citime, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the number of rows in the input file.
    long nrows;
    fits_get_num_rows(infptr, &nrows, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the number of bins for the light curve.
    long nbins=(long)(par.length/par.dt);

    // Allocate memory for the count histogram.
    headas_chat(5, "create empty light curve with %ld bins ...\n",
		nbins);
    counts=(long*)malloc(nbins*sizeof(long));
    CHECK_NULL_BREAK(counts, status, "memory allocation for light curve failed");

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
      fits_read_col(infptr, TDOUBLE, citime, ii+1, 1, 1, 
		    &dnull, &time, &anynul, &status);
      CHECK_STATUS_BREAK(status);
      
      // Determine the respective bin in the light curve.
      long bin=(long)(time/par.dt);

      // If the event exceeds the end of the light curve, simply neglect it.
      if (bin>=nbins) continue;
      
      // Add the event to the light curve.
      assert(bin>=0);
      counts[bin]++;
    }
    CHECK_STATUS_BREAK(status);
    // END of loop over all events in the input file.

    // Store the light curve in the output file.
    headas_chat(3, "store light curve ...\n");

    // Create a new FITS-file (remove existing one before):
    remove(par.LightCurve);
    char buffer[MAXFILENAME];
    sprintf(buffer, "%s(%s%s)", par.LightCurve, SIXT_DATA_PATH, 
	    "/templates/genlc.tpl");
    fits_create_file(&outfptr, buffer, &status);
    CHECK_STATUS_BREAK(status);

    // Move to the HDU containing the binary table.
    int hdutype;
    fits_movabs_hdu(outfptr, 2, &hdutype, &status);
    CHECK_STATUS_BREAK(status);

    // Get column numbers.
    int cotime, ccounts;
    fits_get_colnum(outfptr, CASEINSEN, "TIME", &cotime, &status);
    fits_get_colnum(outfptr, CASEINSEN, "COUNTS", &ccounts, &status);
    CHECK_STATUS_BREAK(status);

    // Write header keywords.
    fits_update_key(outfptr, TSTRING, "TIMEUNIT", "s", 
		    "time unit", &status);
    fits_update_key(outfptr, TDOUBLE, "TIMERES", &par.dt, 
		    "time resolution", &status);
    CHECK_STATUS_BREAK(status);

    // Write the data into the table.
    for (ii=0; ii<nbins; ii++) {
      // Convert the count histogram to a light curve with
      // time and rate entries.
      double dbuffer=(ii+1)*par.dt;
      fits_write_col(outfptr, TDOUBLE, cotime, ii+1, 1, 1, 
		     &dbuffer, &status);
      CHECK_STATUS_BREAK(status);
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

  if (status == EXIT_SUCCESS) headas_chat(3, "finished successfully!\n\n");
  return(status);
}


int genlc_getpar(struct Parameters* par)
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

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  return(status);
}



