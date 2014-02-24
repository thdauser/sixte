/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
*/

#include "badpixmap.h"


BadPixMap* newBadPixMap(int* const status)
{
  headas_chat(3, "initialize empty BadPixMap object ...\n");

  BadPixMap* map=(BadPixMap*)malloc(sizeof(BadPixMap));
  if (NULL==map) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for BadPixMap failed");
    return(map);
  }

  // Initialize all pointers with NULL.
  map->pixels   =NULL;
  map->anybadpix=NULL;

  // Set the initial values.
  map->xwidth=0;
  map->ywidth=0;

  return(map);
}


BadPixMap* loadBadPixMap(const char* const filename, 
			 int* const status)
{
  // Get a new and empty data structure from the constructor.
  BadPixMap* map=NULL;

  fitsfile* fptr=NULL;
  float* input_buffer=NULL; // Input buffer for the bad pixel map image.

  do { // Error handling loop.
    map=newBadPixMap(status);
    CHECK_STATUS_BREAK(*status);

    // Load the data from the FITS file. 
    // Open the FITS file.
    fits_open_image(&fptr, filename, READONLY, status);
    CHECK_STATUS_BREAK(*status);
      
    // Determine the image dimensions.
    long naxes[2];
    fits_get_img_size(fptr, 2, naxes, status);
    CHECK_STATUS_BREAK(*status);
    map->xwidth=(int)naxes[0];
    map->ywidth=(int)naxes[1];

    // Allocate memory for the input buffer.
    input_buffer=(float*)malloc(map->xwidth*map->ywidth*sizeof(float));
    if (NULL==input_buffer) {
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for BadPixMap failed");
      break;
    }
    
    // Read the FITS image.
    int anynul=0;
    int null_value=0.;
    long fpixel[2]={1, 1}; // lower left corner
    //              |--|--> FITS coordinates start at (1,1)
    long lpixel[2]={map->xwidth, map->ywidth}; // upper right corner
    long inc[2]={1, 1};
    fits_read_subset(fptr, TFLOAT, fpixel, lpixel, inc, &null_value,
		     input_buffer, &anynul, status);
    CHECK_STATUS_BREAK(*status);
    
    // Allocate memory for the pixel array.
    map->pixels=(float**)malloc(map->xwidth*sizeof(float*));
    if (NULL==map->pixels) {
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for BadPixMap failed");
      break;
    }
    int ii;
    for (ii=0; ii<map->xwidth; ii++) {
      map->pixels[ii]=(float*)malloc(map->ywidth*sizeof(float));
      if (NULL==map->pixels[ii]) {
	*status=EXIT_FAILURE;
	SIXT_ERROR("memory allocation for BadPixMap failed");
	break;
      }
    }
    CHECK_STATUS_BREAK(*status);
    // Allocate memory for the column flags.
    map->anybadpix=(int*)malloc(map->xwidth*sizeof(int));
    if (NULL==map->anybadpix) {
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for BadPixMap failed");
      break;
    }
    for (ii=0; ii<map->xwidth; ii++) {
      map->anybadpix[ii]=0;
    }

    // Copy the input buffer to the BadPixMap.
    for(ii=0; ii<map->xwidth; ii++) {
      int jj;
      for(jj=0; jj<map->ywidth; jj++) {
        map->pixels[ii][jj]=input_buffer[ii+ map->xwidth*jj];
        // Check if the pixel is a bad one.
        if (map->pixels[ii][jj]!=0.) {
	  // Set the flag that this column of the image contains
	  // at least one bad pixel.
	  map->anybadpix[ii]=1;

	  // Print some information.
	  if (map->pixels[ii][jj]>0.) {
	    headas_chat(5, " hot pixel at (%d,%d)\n", ii, jj);
	  } else if (map->pixels[ii][jj]<0.) {
	    headas_chat(5, " cold pixel at (%d,%d)\n", ii, jj);
	  }
        }
      }
    }
  } while (0); // END of 1st error handling loop.

  // Close the FITS file.
  if (NULL!=fptr) {
    fits_close_file(fptr, status);
  }

  // Clean up.
  if (NULL!=input_buffer) {
    free(input_buffer);
  }

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
