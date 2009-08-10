#include "hexagonalpixels.h"

int initHexagonalPixels(HexagonalPixels* hp, struct HexagonalPixelsParameters* parameters) 
{
  int status = EXIT_SUCCESS;

  // Set the length of the pixel array:
  hp->npixels = parameters->npixels;
  // Set the pixel dimensions:
  hp->h = parameters->pixelwidth/2.;
  hp->a = parameters->pixelwidth/sqrt(3.);

  // Get the memory for the pixels:
  hp->array = (HexagonalPixel*)malloc(hp->npixels*sizeof(HexagonalPixel));
  if (NULL==hp->array) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: memory allocation for pixel array failed!\n", status);
    return(status);
  }

  // Clear the pixels.
  //clearHexagonalPixels(sp);

  return(status);
}



void cleanupHexagonalPixels(HexagonalPixels* hp) 
{
  if (NULL!=hp->array) {
    free(hp->array);
    hp->array=NULL;
  }
}



inline void clearHexagonalPixels(HexagonalPixels* hp)
{
  int x;
  for (x=0; x<hp->npixels; x++) {
    hp->array[x].charge = 0.;
  }
}



