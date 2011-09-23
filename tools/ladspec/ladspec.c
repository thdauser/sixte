#include "ladspec.h"


////////////////////////////////////
/** Main procedure. */
int ladspec_main() {
  // Program parameters.
  struct Parameters par; 

  // TODO Use variable number of bins and energy range,
  // which can be selected by the user.

  // Number of bins in the spectrum
  const long NBINS=1024;


  // Input event list file.
  LADEventListFile* elf=NULL;

  // Output spectrum.
  long* spec=NULL;
  fitsfile* fptr=NULL;

  // Error status.
  int status=EXIT_SUCCESS;   


  // Register HEATOOL:
  set_toolname("ladspec");
  set_toolversion("0.01");


  do {  // Beginning of the ERROR handling loop.

    // --- Initialization ---

    // Read the program parameters using PIL library.
    status=ladspec_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Set the input pattern file.
    elf=openLADEventListFile(par.EventList, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Allocate memory for the output spectrum.
    spec=(long*)malloc(NBINS*sizeof(long));
    CHECK_NULL_BREAK(spec, status, "memory allocation for spectrum failed");

    // Initialize the spectrum with 0.
    long ii; 
    for (ii=0; ii<NBINS; ii++) {
      spec[ii]=0;
    }

    // --- END of Initialization ---


    // --- Begin Spectrum Binning ---
    headas_chat(5, "create spectrum ...\n");

    // LOOP over all events in the FITS table.
    for (ii=0; ii<elf->nrows; ii++) {
      
      // Read the next event from the file.
      LADEvent event;
      getLADEventFromFile(elf, ii+1, &event, &status);
      CHECK_STATUS_BREAK(status);
      
    }
    CHECK_STATUS_BREAK(status);
    // END of loop over all events in the input file.

    // Store the spectrum in the output file.
    headas_chat(3, "store spectrum ...\n");

    // Create a new FITS-file (remove existing one before):
    remove(par.Spectrum);
    fits_create_file(&fptr, par.Spectrum, &status);
    CHECK_STATUS_BREAK(status);

    // TODO

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Close the files.
  freeLADEventListFile(&elf, &status);
  if (NULL!=fptr) fits_close_file(fptr, &status);

  // Release memory.
  if (NULL!=spec) free(spec);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully!\n\n");
  return(status);
}


int ladspec_getpar(struct Parameters* par)
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

  status=ape_trad_query_file_name("Spectrum", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the output spectrum file");
    return(status);
  } 
  strcpy(par->Spectrum, sbuffer);
  free(sbuffer);

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  return(status);
}



