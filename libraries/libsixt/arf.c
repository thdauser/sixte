#include "arf.h"


struct ARF* loadARF(const char* const filename, int* const status) 
{
  fitsfile* fptr=NULL;

  struct ARF* arf = (struct ARF*)malloc(sizeof(struct ARF));
  if (NULL==arf) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: could not allocate memory for ARF!\n", *status);
    return(arf);
  }

  // Load the RMF from the FITS file using the HEAdas RMF access routines
  // (part of libhdsp).
  fits_open_file(&fptr, filename, READONLY, status);
  if (*status!=EXIT_SUCCESS) return(arf);
  
  // Read the 'SPECRESP MATRIX' or 'MATRIX' extension:
  if ((*status=ReadARF(fptr, 0, arf))!=EXIT_SUCCESS) return(arf);

  // Print some information:
  headas_chat(5, "ARF loaded with %ld energy bins\n",
	      arf->NumberEnergyBins);

  // Close the open FITS file.
  fits_close_file(fptr, status);
  
  return(arf);
}


void freeARF(struct ARF* const arf) 
{
  if (NULL!=arf) {
    free(arf);
  }
}

