#include "arcpixels.h"


int initArcPixels(ArcPixels* ap, struct ArcPixelsParameters* app)
{
  int ring;  // Counter for the individual detector rings.
  int status = EXIT_SUCCESS;

  // Set the array/detector dimensions:
  ap->nrings = app->nrings;
  ap->npixels = (int*)malloc(ap->nrings*sizeof(int));
  if (NULL==ap->npixels) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: memory allocation for ArcPixels failed!\n", status);
    return(status);
  }
  ap->radius = (double*)malloc(ap->nrings*sizeof(double));
  if (NULL==ap->radius) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: memory allocation for ArcPixels failed!\n", status);
    return(status);
  }
  ap->offset_angle = (double*)malloc(ap->nrings*sizeof(double));
  if (NULL==ap->offset_angle) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: memory allocation for ArcPixels failed!\n", status);
    return(status);
  }
  for (ring=0; ring<ap->nrings; ring++) {
    ap->npixels[ring] = app->npixels[ring];
    ap->radius[ring]  = app->radius[ring];
    ap->offset_angle[ring] = app->offset_angle[ring];
  }

  // Get the memory for the pixels:
  ap->array = (ArcPixel**)malloc(ap->nrings*sizeof(ArcPixel*));
  if (NULL==ap->array) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: memory allocation for pixel array failed!\n", status);
    return(status);
  }
  for(ring=0; ring<ap->nrings; ring++) {
    ap->array[ring] = (ArcPixel*)malloc(ap->npixels[ring]*sizeof(ArcPixel));
    if (NULL==ap->array[ring]) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for pixel array failed!\n", status);
      return(status);
    }
  }

  // Clear the pixels.
  clearArcPixels(ap);

  return(status);
}



inline void clearArcPixels(ArcPixels* ap) 
{
  int ring, pixel;
  for (ring=0; ring<ap->nrings; ring++) {
    for (pixel=0; pixel<ap->npixels[ring]; pixel++) {
      ap->array[ring][pixel].charge = 0.;
      ap->array[ring][pixel].last_impact = -1.e20;
    }
  }
}



void cleanupArcPixels(ArcPixels* ap) 
{
  // Free the array/detector dimensions.
  if (NULL!=ap->npixels) {
    free(ap->npixels);
    ap->npixels=NULL;
  }
  if (NULL!=ap->radius) {
    free(ap->radius);
    ap->radius=NULL;
  }
  if (NULL!=ap->offset_angle) {
    free(ap->offset_angle);
    ap->offset_angle=NULL;
  }

  // Free the pixel array.
  if (NULL!=ap->array) {
    int ring;
    for (ring=0; ring<ap->nrings; ring++) {
      if (NULL!=ap->array[ring]) {
	free(ap->array[ring]);
      }
    }
    free(ap->array);
    ap->array=NULL;
  }
}



void getArcPixel(ArcPixels* ap, struct Point2d position, int* ring, int* number)
{
  // Search for the corresponding pixel ring.
  double radius = sqrt(pow(position.x,2.)+pow(position.y,2.));
  for(*ring=0; *ring<ap->nrings; (*ring)++) {
    if (radius < ap->radius[*ring]) break;
  }
  if (*ring==ap->nrings) {
    *ring=INVALID_PIXEL;
    *number=INVALID_PIXEL;
    return;
  }

  // Search for the pixel index within this ring.
  double angle = atan2(position.y, position.x); // Angle is within [-pi:pi].
  angle -= ap->offset_angle[*ring];  // Subtract the offset angle.
  if (angle < 0.) angle+=2.*M_PI;   // Now the angle is within [0:2pi]
  *number = (int)(angle/(2.*M_PI/ap->npixels[*ring]));
}



int getArcPixelIndex(ArcPixels* ap, int ring, int number)
{
  int count, index=0;
  for(count=0; count<ring; count++) {
    index += ap->npixels[count];
  }
  return(index+number);
}

