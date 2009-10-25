#include "codedmask.h"


CodedMask* getCodedMask(int* status)
{
  CodedMask* mask=(CodedMask*)malloc(sizeof(CodedMask));
  if (NULL==mask) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for Coded Mask!\n",
		   *status);
    return(mask);
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



CodedMask* getCodedMaskFromFile(char* filename, int* status)
{
  // Obtain a new (empty) CodedMask object using the basic constructor.
  CodedMask* mask = getCodedMask(status);
  if (EXIT_SUCCESS!=*status) return(mask);

  fitsfile* fptr=NULL;
  int* input_buffer=NULL; // Input buffer for FITS image.

  do { // Beginning of Error Handling loop. 

    // Open the FITS file containing the coded mask.
    headas_chat(5, "open Coded Mask FITS file '%s' ...\n", filename);
    if (fits_open_image(&fptr, filename, READONLY, status)) break;

    // Determine the width of the image.
    long naxes[2];
    if (fits_get_img_size(fptr, 2, naxes, status)) break;
    assert(naxes[0]==naxes[1]);
    mask->naxis1 = (int)naxes[0];
    mask->naxis2 = (int)naxes[1];
    
    // Determine the width of one image pixel.
    char comment[MAXMSG]; // buffer
    if (fits_read_key(fptr, TDOUBLE, "CDELT1", &mask->cdelt1, comment, status))
      break;
    if (fits_read_key(fptr, TDOUBLE, "CDELT2", &mask->cdelt2, comment, status))
      break;
    // Determine the WCS keywords of the image.
    if (fits_read_key(fptr, TDOUBLE, "CRPIX1", &mask->crpix1, comment, status))
      break;
    if (fits_read_key(fptr, TDOUBLE, "CRPIX2", &mask->crpix2, comment, status))
      break;
    if (fits_read_key(fptr, TDOUBLE, "CRVAL1", &mask->crval1, comment, status))
      break;
    if (fits_read_key(fptr, TDOUBLE, "CRVAL2", &mask->crval2, comment, status))
      break;

    // Allocate memory for the pixels of the image:
    mask->map = (int**)malloc(mask->naxis1*sizeof(int*));
    if (NULL!=mask->map) {
      int count;
      for(count=0; (count<mask->naxis1)&&(EXIT_SUCCESS==*status); count++) {
	mask->map[count] = (int*)malloc(mask->naxis2*sizeof(int));
	if(NULL==mask->map[count]) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: could not allocate memory to store the "
			 "CodedMask!\n", *status);
	  break;
	}
      }
    } else {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store the "
		     "CodedMask!\n", *status);
      break;
    }     
    // Allocate memory for input buffer (1D array):
    input_buffer=(int*)malloc(mask->naxis1*mask->naxis2*sizeof(int));
    if(NULL==input_buffer) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store the "
		     "CodedMask!\n", *status);
      break;
    }
    // END of memory allocation


  } while (0); // END of Error Handling loop.

  // Release memory from image input buffer.
  if (NULL!=input_buffer) free(input_buffer);

  // Close the FITS file.
  if (NULL!=fptr) fits_close_file(fptr, status);

  return(mask);
}

