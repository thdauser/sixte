#include "gentel.h"


////////////////////////////////////////////////////////////////////
// Program Code
////////////////////////////////////////////////////////////////////


GenTel* newGenTel(int* const status) 
{
  // Allocate memory.
  GenTel* tel=(GenTel*)malloc(sizeof(GenTel));
  if (NULL==tel) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for GenTel failed");
    return(tel);
  }

  // Initialize all pointers with NULL.
  tel->psf       =NULL;
  tel->vignetting=NULL;
  tel->coded_mask=NULL;
  tel->arf       =NULL;
  tel->arf_filename=NULL;

  // Set initial values.
  tel->focal_length=0.;
  tel->fov_diameter=0.;

  return(tel);
}


void destroyGenTel(GenTel** const tel)
{
  if (NULL!=*tel) {
    if (NULL!=(*tel)->arf) {
      free((*tel)->arf);
    }
    if (NULL!=(*tel)->arf_filename) {
      free((*tel)->arf_filename);
    }
    destroyPSF(&(*tel)->psf);
    destroyVignetting(&(*tel)->vignetting);
    destroyCodedMask(&(*tel)->coded_mask);

    free(*tel);
    *tel=NULL;
  }
}


