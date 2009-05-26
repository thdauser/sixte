/** Data structure containing the mirror vignetting function. */
typedef struct {
  float* data;
} Vignetting;



/** Constructor of the Vignetting data structure. 
 * Loads the vignetting function from a given FITS file. 
 * The format of the FITS file is defined by 
 * OGIP Memo CAL/GEN/92-021. */
Vignetting* get_Vignetting(char* filename, int* status) {
  Vignetting* vignetting=NULL;
  fitsfile* fptr=NULL;

  do {

    // Open the FITS file for reading the vignetting function.
    if (fits_open_table(&fptr, filename, READONLY, status)) break;

  } while(0); // END of Error handling loop

  // --- Clean up ---

  if (NULL!=fptr) fits_close_file(fptr, status);

  if (EXIT_SUCCESS!=*status) vignetting=NULL;
  return(vignetting);
}


