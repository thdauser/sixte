#include "ladlc.h"


////////////////////////////////////
/** Main procedure. */
int ladlc_main() {
  // Program parameters.
  struct Parameters par; 

  // Input event list file.
  LADEventListFile* elf=NULL;

  // Output light curve.
  long* counts=NULL;
  fitsfile* fptr=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL:
  set_toolname("ladlc");
  set_toolversion("0.01");


  do {  // Beginning of the ERROR handling loop.

    // --- Initialization ---

    // Read the program parameters using PIL library.
    status=ladlc_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Set the input pattern file.
    elf=openLADEventListFile(par.EventList, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the number of bins for the light curve.
    long nbins=(long)(par.length/par.dt)+1;

    // Allocate memory for the output light curve.
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
    headas_chat(5, "calculate light curve ...\n");

    // LOOP over all events in the FITS table.
    for (ii=0; ii<elf->nrows; ii++) {
      
      // Read the next event from the file.
      LADEvent event;
      getLADEventFromFile(elf, ii+1, &event, &status);
      CHECK_STATUS_BREAK(status);
      
      // Determine the respective bin in the light curve.
      long bin=(long)(event.time/par.dt);
      
      // Add the event to the light curve.
      assert(bin>=0);
      assert(bin<nbins);      
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
	    "/templates/ladlc.tpl");
    fits_create_file(&fptr, buffer, &status);
    CHECK_STATUS_BREAK(status);

    // Move to the HDU containing the binary table.
    int hdutype;
    fits_movabs_hdu(fptr, 2, &hdutype, &status);
    CHECK_STATUS_BREAK(status);

    // Get column numbers.
    int ccounts;
    fits_get_colnum(fptr, CASEINSEN, "COUNTS", &ccounts, &status);
    CHECK_STATUS_BREAK(status);

    // Write header keywords.
    fits_update_key(fptr, TDOUBLE, "TIMERES", &par.dt, 
		    "time resolution", &status);
    CHECK_STATUS_BREAK(status);

    // Write the data into the table.
    fits_write_col(fptr, TLONG, ccounts, 1, 1, nbins, counts, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Close the files.
  freeLADEventListFile(&elf, &status);
  if (NULL!=fptr) fits_close_file(fptr, &status);

  // Release memory.
  if (NULL!=counts) free(counts);

  if (status == EXIT_SUCCESS) headas_chat(3, "finished successfully!\n\n");
  return(status);
}


int ladlc_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status = EXIT_SUCCESS; 

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



