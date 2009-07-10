#include "squarepixels.h"


int initSquarePixels(SquarePixels* sp, struct SquarePixelsParameters* parameters) 
{
  int status = EXIT_SUCCESS;

  // Set the array dimensions:
  sp->xwidth = parameters->xwidth;
  sp->ywidth = parameters->ywidth;
  sp->xoffset = parameters->xwidth/2;
  sp->yoffset = parameters->ywidth/2;
  sp->xpixelwidth = parameters->xpixelwidth;
  sp->ypixelwidth = parameters->ypixelwidth;

  // Get the memory for the pixels:
  sp->array = (SquarePixel**)malloc(sp->xwidth*sizeof(SquarePixel*));
  if (NULL==sp->array) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: memory allocation for pixel array failed!\n", status);
    return(status);
  } else {
    int count;
    for(count=0; count<sp->xwidth; count++) {
      sp->array[count] = (SquarePixel*)malloc(sp->ywidth*sizeof(SquarePixel));
      if (NULL==sp->array[count]) {
	status = EXIT_FAILURE;
	HD_ERROR_THROW("Error: memory allocation for pixel array failed!\n", status);
	return(status);
      }
    }
  }

  return(status);
}



