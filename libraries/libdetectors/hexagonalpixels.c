#include "hexagonalpixels.h"

int initHexagonalPixels(HexagonalPixels* hp, struct HexagonalPixelsParameters* parameters) 
{
  int status = EXIT_SUCCESS;

  // Set the length of the pixel array:
  hp->npixels = parameters->npixels;
  assert(HTRS_N_PIXELS==hp->npixels);

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
  clearHexagonalPixels(hp);



  // Set up some auxiliary data that is required to determine the affected pixels
  // for a cartesian photon impact position.
  // Determine the coordinates of the pixel centers:
  struct Point2d centers[HTRS_N_PIXELS];

  centers[0].x  =  0.; // Pixel number 1
  centers[0].y  =  0.; 
  //
  centers[1].x  = -1.5 * hp->a;
  centers[1].y  = -1.  * hp->h;

  centers[2].x  =  0.;
  centers[2].y  = -2.  * hp->h;

  centers[3].x  =  1.5 * hp->a;
  centers[3].y  = -1.  * hp->h;

  centers[4].x  =  1.5 * hp->a;
  centers[4].y  =  1.  * hp->h;

  centers[5].x  =  0.;
  centers[5].y  =  2.  * hp->h;

  centers[6].x  = -1.5 * hp->a;
  centers[6].y  =  1.  * hp->h;
  //
  centers[7].x  = -3.  * hp->a;
  centers[7].y  =  0.;

  centers[8].x  = -3.  * hp->a;
  centers[8].y  = -2.  * hp->h;

  centers[9].x  = -1.5 * hp->a;
  centers[9].y  = -3.  * hp->h;

  centers[10].x =  0.;
  centers[10].y = -4.  * hp->h;

  centers[11].x =  1.5 * hp->a;
  centers[11].y = -3.  * hp->h;

  centers[12].x =  3.  * hp->a;
  centers[12].y = -2.  * hp->h;

  centers[13].x =  3.  * hp->a;
  centers[13].y =  0.;

  centers[14].x =  3.  * hp->a;
  centers[14].y =  2.  * hp->h;

  centers[15].x =  1.5 * hp->a;
  centers[15].y =  3.  * hp->h;

  centers[16].x =  0.;
  centers[16].y =  4.  * hp->h;

  centers[17].x = -1.5 * hp->a;
  centers[17].y =  3.  * hp->h;

  centers[18].x = -3.  * hp->a;
  centers[18].y =  2.  * hp->h;
  //
  centers[19].x = -4.5 * hp->a;
  centers[19].y =  1.  * hp->h;

  centers[20].x = -4.5 * hp->a;
  centers[20].y = -1.  * hp->h;

  centers[21].x = -4.5 * hp->a;
  centers[21].y = -3.  * hp->h;

  centers[22].x = -3.  * hp->a;
  centers[22].y = -4.  * hp->h;

  centers[23].x = -1.5 * hp->a;
  centers[23].y = -5.  * hp->h;

  centers[24].x =  0.;
  centers[24].y = -6.  * hp->h;

  centers[25].x =  1.5 * hp->a;
  centers[25].y = -5.  * hp->h;

  centers[26].x =  3.  * hp->a;
  centers[26].y = -4.  * hp->h;

  centers[27].x =  4.5 * hp->a;
  centers[27].y = -3.  * hp->h;

  centers[28].x =  4.5 * hp->a;
  centers[28].y = -1.  * hp->h;

  centers[29].x =  4.5 * hp->a;
  centers[29].y =  1.  * hp->h;

  centers[30].x =  4.5 * hp->a;
  centers[30].y =  3.  * hp->h;

  centers[31].x =  3.  * hp->a;
  centers[31].y =  4.  * hp->h;

  centers[32].x =  1.5 * hp->a;
  centers[32].y =  5.  * hp->h;

  centers[33].x =  0.;
  centers[33].y =  6.  * hp->h;

  centers[34].x = -1.5 * hp->a;
  centers[34].y =  5.  * hp->h;

  centers[35].x = -3.  * hp->a;
  centers[35].y =  4.  * hp->h;

  centers[36].x = -4.5 * hp->a;
  centers[36].y =  3.  * hp->h;



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



