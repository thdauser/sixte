#include "codedmask.h"


CodedMask* getCodedMask(int* status)
{
  CodedMask* mask=(CodedMask*)malloc(sizeof(CodedMask));
  if (NULL==mask) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for Coded Mask!\n",
		   *status);
    return(NULL);
  }

  // Set initial default values.
  mask->map=NULL;
  mask->transparent_pixels=NULL;
  mask->n_transparent_pixels=0;
  mask->opacity=0.;
  mask->naxis1=0;
  mask->naxis2=0;
  mask->cdelt1=0.;
  mask->cdelt2=0.;
  mask->crpix1=0.;
  mask->crpix2=0.;
  mask->crval1=0.;
  mask->crval2=0.;

  return(mask);
}



void freeCodedMask(CodedMask* mask)
{
  if (NULL!=mask) {
    if (NULL!=mask->map) {
      int count;
      for(count=0; count<mask->naxis1; count++) {
	if (NULL!=mask->map[count]) {
	  free(mask->map[count]);
	}
      }
      free(mask->map);
    }
    if (NULL!=mask->transparent_pixels) {
      int count;
      for(count=0; count<mask->n_transparent_pixels; count++) {
	if (NULL!=mask->transparent_pixels[count]) {
	  free(mask->transparent_pixels[count]);
	}
      }
      free(mask->transparent_pixels);
    }
    free(mask);
  }
}

