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
  mask->transparency=0.;
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
  int x, y;

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
      for(x=0; x<mask->naxis1; x++) {
	mask->map[x] = (int*)malloc(mask->naxis2*sizeof(int));
	if(NULL==mask->map[x]) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: could not allocate memory to store the "
			 "CodedMask!\n", *status);
	  break;
	}
      }
      if (EXIT_SUCCESS!=*status) break;
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

    
    // READ the FITS image:
    int anynul=0;
    int null_value=0.;
    long fpixel[2] = {1, 1}; // lower left corner
    //                |--|--> FITS coordinates start at (1,1)
    long lpixel[2] = {mask->naxis1, mask->naxis2}; // upper right corner
    long inc[2] = {1, 1};
    if (fits_read_subset(fptr, TINT, fpixel, lpixel, inc, &null_value, 
			 input_buffer, &anynul, status)) break;

    // Copy the image from the input buffer to the map array.
    mask->n_transparent_pixels=0;
    for(x=0; x<mask->naxis1; x++) {
      for(y=0; y<mask->naxis2; y++) {
	mask->map[x][y] = input_buffer[x+ mask->naxis2*y];
	
	// Check if the pixel is transparent or opaque.
	if (TRANSPARENT==mask->map[x][y]) {
	  mask->n_transparent_pixels++;
	}
      }
    }


    // Determine the masks transparency.
    mask->transparency = 
      ((double)mask->n_transparent_pixels)/((double)(mask->naxis1*mask->naxis2));
    
    // Allocate memory for the list of transparent pixels.
    mask->transparent_pixels = (int**)malloc(mask->n_transparent_pixels*sizeof(int*));
    if (NULL==mask->transparent_pixels) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store the "
		     "CodedMask!\n", *status);
      break;
    }
    for(x=0; x<mask->n_transparent_pixels; x++) {
      mask->transparent_pixels[x] = (int*)malloc(2*sizeof(int));
      if (NULL==mask->transparent_pixels[x]) {
	*status=EXIT_FAILURE;
	HD_ERROR_THROW("Error: could not allocate memory to store the "
		       "CodedMask!\n", *status);
	break;
      }
    }
    if (EXIT_SUCCESS!=*status) break;

    int count=0;
    for (x=0; x<mask->naxis1; x++) {
      for (y=0; y<mask->naxis2; y++) {
	if (TRANSPARENT==mask->map[x][y]) {
	  mask->transparent_pixels[count][0] = x;
	  mask->transparent_pixels[count][1] = y;
	  count++;
	}
      }
    }  

  } while (0); // END of Error Handling loop.

  // Release memory from image input buffer.
  if (NULL!=input_buffer) free(input_buffer);

  // Close the FITS file.
  if (NULL!=fptr) fits_close_file(fptr, status);

  return(mask);
}



int getCodedMaskImpactPos(struct Point2d* position, Photon* photon, CodedMask* mask, 
			  struct Telescope* telescope)
{

  // Get a random number.
  double rand = get_random_number();

  if (rand>mask->transparency) {
    // The photon is absorbed by an opaque pixel.
    return(0);
  }

  // The photon passes through a tranparent pixel. So we have to 
  // determine its impact position on the detector plane using
  // geometrical considerations.
  // First determine the pixel where the photon passes.
  int pixel = (int)(rand/mask->transparency * mask->n_transparent_pixels);
  // Now pixel points to an arbitrary pixel in the list of transparent pixels.

  // Determine the direction of origin of the photon.
  Vector photon_direction = unit_vector(photon->ra, photon->dec);
  // Determine the components of the photon direction vector with respect to the
  // detector coordinate axes nx and ny.
  double scpx = scalar_product(&photon_direction, &telescope->nx);
  double scpy = scalar_product(&photon_direction, &telescope->ny);
  // Determine the component of the photon direction within the detector plane.
  double radius = sqrt(pow(scpx,2.)+pow(scpy,2.));
  // Determine the azimuthal angle of the photon within the detector plane
  // with respect to the detector nx axis.
  double alpha = atan2(scpy, scpx);

  // Determine the impact position in the mask plane.
  position->x = ((double)(mask->transparent_pixels[pixel][0])
		 -mask->crpix1+0.5+get_random_number())*mask->cdelt1 + mask->crval1;
  position->y = ((double)(mask->transparent_pixels[pixel][1])
		 -mask->crpix2+0.5+get_random_number())*mask->cdelt2 + mask->crval2;

  // Shift the position to obtain the position in the detector plane.
  // The shift is necessary because of the photon's off-axis position.
  position->x -= cos(alpha) * radius * telescope->focal_length;
  position->y -= sin(alpha) * radius * telescope->focal_length;

  return(1);
}


