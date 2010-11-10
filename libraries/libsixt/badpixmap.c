#include "badpixmap.h"


BadPixMap* newBadPixMap(int* const status)
{
  headas_chat(5, "initialize empty BadPixMap object ...\n");

  BadPixMap* map = (BadPixMap*)malloc(sizeof(BadPixMap));
  if (NULL==map) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for BadPixMap failed!\n", *status);
    return(map);
  }

  // Initialize all pointers with NULL.
  map->pixels    = NULL;
  map->anybadpix = NULL;

  // Set the initial values.
  map->xwidth = 0;
  map->ywidth = 0;

  return(map);
}


BadPixMap* loadBadPixMap(const char* const filename, 
			 int* const status)
{
  // Get a new and empty data structure from the constructor.
  BadPixMap* map = newBadPixMap(status);
  if (EXIT_SUCCESS!=*status) return(map);

  // Load the data from the FITS file. 
  float* input_buffer=NULL; // Input buffer for the bad pixel map image.

  do { // 1st ror handling loop.

    fitsfile* fptr=NULL;
    do { // 2nd error handling loop.
      
      // Open the FITS file.
      if (fits_open_image(&fptr, filename, READONLY, status)) break;
      
      // Determine the image dimensions.
      long naxes[2];
      if (fits_get_img_size(fptr, 2, naxes, status)) break;
      map->xwidth = (int)naxes[0];
      map->ywidth = (int)naxes[1];

      // Allocate memory for the input buffer.
      input_buffer=(float*)malloc(map->xwidth*map->ywidth*sizeof(float));
      if(NULL==input_buffer) {
	*status=EXIT_FAILURE;
	HD_ERROR_THROW("Error: could not allocate memory to store the "
		       "BadPixMap!\n", *status);
	break;
      }

      // Read the FITS image.
      int anynul=0;
      int null_value=0.;
      long fpixel[2] = {1, 1}; // lower left corner
      //                |--|--> FITS coordinates start at (1,1)
      long lpixel[2] = {map->xwidth, map->ywidth}; // upper right corner
      long inc[2] = {1, 1};
      if (fits_read_subset(fptr, TFLOAT, fpixel, lpixel, inc, &null_value,
			   input_buffer, &anynul, status)) break;
      
    } while (0); // END of 2nd error handling loop.
    
    // Close the FITS file.
    fits_close_file(fptr, status);
    if (EXIT_SUCCESS!=*status) break;
    
    // Allocate memory for the pixel array.
    map->pixels = (float**)malloc(map->xwidth*sizeof(float*));
    if (NULL==map->pixels) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: Memory allocation for BadPixMap failed!\n", *status);
      return(map);
    }
    int ii;
    for (ii=0; ii<map->xwidth; ii++) {
      map->pixels[ii] = (float*)malloc(map->ywidth*sizeof(float));
      if (NULL==map->pixels[ii]) {
	*status=EXIT_FAILURE;
	HD_ERROR_THROW("Error: Memory allocation for BadPixMap failed!\n", *status);
	return(map);
      }
    }
    // Allocate memory for the column flags.
    map->anybadpix = (int*)malloc(map->xwidth*sizeof(int));
    if (NULL==map->anybadpix) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: Memory allocation for BadPixMap failed!\n", *status);
      return(map);
    }
    for (ii=0; ii<map->xwidth; ii++) {
      map->anybadpix[ii] = 0;
    }

    // Copy the input buffer to the BadPixMap.
    for(ii=0; ii<map->xwidth; ii++) {
      int jj;
      for(jj=0; jj<map->ywidth; jj++) {
        map->pixels[ii][jj] = input_buffer[ii+ map->xwidth*jj];
        // Check if the pixel is a bad one.
        if (map->pixels[ii][jj]!=0.) {
	  // Set the flag that this column of the image contains
	  // at least one bad pixel.
	  map->anybadpix[ii] = 1;
        }

	// Print some information.
	if (map->pixels[ii][jj]>0.) {
	  headas_chat(4, " hot pixel at (%d,%d)\n", ii, jj);
	} else if (map->pixels[ii][jj]<0.) {
	  headas_chat(4, " cold pixel at (%d,%d)\n", ii, jj);
	}
      }
    }

  } while (0); // END of 1st error handling loop.

  // Clean up.
  if (NULL!=input_buffer) {
    free(input_buffer);
  }
  if (EXIT_SUCCESS!=*status) return(map);
  // END of loading the data from the FITS file.


  return(map);
}


void destroyBadPixMap(BadPixMap** const map) 
{
  if (NULL!=(*map)) {
    if (NULL!=(*map)->pixels) {
      int ii;
      for (ii=0; ii<(*map)->xwidth; ii++) {
	if (NULL!=(*map)->pixels[ii]) {
	  free((*map)->pixels[ii]);
	}
      }
      free((*map)->pixels);
    }
    if (NULL!=(*map)->anybadpix) {
      free((*map)->anybadpix);
    }
    free(*map);
    *map=NULL;
  }
}



void applyBadPixMap(const BadPixMap* const map, const double timespan,
		    void (*encounter) (void* const data, 
				       const int x, const int y, 
				       const float value),
		    void* const data)
{
  int ii;
  for (ii=0; ii<map->xwidth; ii++) {
    if (1==map->anybadpix[ii]) {
      int jj;
      for (jj=0; jj<map->ywidth; jj++) {
	if (0.!=map->pixels[ii][jj]) {
	  // Call the encounter routine with the value of the bad pixel
	  // multiplied with the timespan weigth.
	  (*encounter)(data, ii, jj, map->pixels[ii][jj]*timespan);
	}
      }
      // END of loop over y-coordinate.
    }
  } 
  // END of loop over x-coordinate.
}


