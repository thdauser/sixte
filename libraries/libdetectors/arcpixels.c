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



int getArcPixelSplits(ArcPixels* ap, GenericDetector* gd, 
		      struct Point2d position, 
		      int* ring, int* number)
{
  // Exact impact position in polar coordinates.
  double radius, angle;

  // Determine the polar coordinates of the impact position.
  getPolarCoordinates(position, &radius, &angle);

  // Determine the main impact position.
  getArcPixelFromPolar(ap, radius, angle, &(ring[0]), &(number[0]));

  // If the impact position is not within a valid pixel, return
  // immediately.
  if (INVALID_PIXEL==ring[0]) return(0);

  // Check if charge cloud size is zero. If this is the case,
  // the event must be a single.
  if (gd->gcc.ccsize < 1.e-20) {
    // Only single events are possible!
    return(1);
  }

  // The event could be a split event (up to double) but also 
  // may be a single.
  
  // Check the 4 next-nearest neighbors:
  // Vary radius:
  if (ring[0] > 0) {
    getArcPixelFromPolar(ap, radius-gd->gcc.ccsize, angle,
			 &(ring[1]), &(number[1]));
    if (ring[1]!=ring[0]) return(2);
  }
  if (ring[0] < ap->nrings-1) {
    getArcPixelFromPolar(ap, radius+gd->gcc.ccsize, angle,
			 &(ring[1]), &(number[1]));
    if (ring[1]!=ring[0]) return(2);
  }
  // Vary polar angle:
  if (ring[0] > 0) {
    double delta = asin(gd->gcc.ccsize*0.5/radius);
    // Bigger polar angle:
    getArcPixelFromPolar(ap, radius, angle+2.*delta,
			 &(ring[1]), &(number[1]));
    if (number[1]!=number[0]) return(2);
    // Smaller polar angle:
    getArcPixelFromPolar(ap, radius, angle-2.*delta,
			 &(ring[1]), &(number[1]));
    if (number[1]!=number[0]) return(2);    
  }
  
  // There is no charge splitting. The event is a single.
  return(1);
}



void getArcPixelFromPolar(ArcPixels* ap, 
			  double radius, double angle,
			  int* ring, int* number)
{
  // Search for the corresponding pixel ring.
  for(*ring=0; *ring<ap->nrings; (*ring)++) {
    if (radius < ap->radius[*ring]) break;
  }
  // Check if the position exceeds the outer detector ring.
  if (*ring==ap->nrings) {
    *ring=INVALID_PIXEL;
    *number=INVALID_PIXEL;
    return;
  }

  // Search for the pixel index within this ring.
  // Subtract the offset angle of the particular pixel ring.
  angle -= ap->offset_angle[*ring];  
  // Make shure that the angle is within [0:2pi].
  while (angle < 0.)       angle+=2.*M_PI;
  while (angle >= 2.*M_PI) angle-=2.*M_PI;
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



void getPolarCoordinates(struct Point2d position, 
			 double* radius, double* angle) 
{
  // Determine the radial distance from the origin.
  *radius = sqrt(pow(position.x,2.)+pow(position.y,2.));

  // Determine the azimuthal angle.
  *angle = atan2(position.y, position.x); // Angle is within [-pi:pi].
  if (*angle < 0.) *angle+=2.*M_PI;   // Now the angle is within [0:2pi].
}

