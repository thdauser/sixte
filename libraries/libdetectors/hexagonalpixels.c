#include "hexagonalpixels.h"


/** Determine the line index l of the line (from the band of parallel lines),
 * where   x > y*m + t_l   and   x < y*m + t_(l+1)
 * with x and y the cartesian coordinates of the photon impact position.
 * The zero line k = 0 with t_k = 0 goes through the center of the detector. */
static inline void getHexagonalPixelLineIndex(struct Point2d position, 
					      double m, double dt, int* l)
{
  // Distance from the zero line of this band of parallel lines.
  double dx = position.x - m*position.y;

  // Determine the resulting line index l.
  // (Avoid integer casting errors by shifting the line.)
  *l = (int)(dx/dt + 100.) - 100;
}



/** Determine the 3 lines that define the sub-triangle of the hexagonal pixel
 * where the impact position is located. */
static inline void getHexagonalPixelLineIndices(HexagonalPixels* hp, struct Point2d position,
						int* l0, int* l1, int* l2)
{
  getHexagonalPixelLineIndex(position,        0.,    hp->h, l0);
  getHexagonalPixelLineIndex(position,  sqrt(3.), 2.*hp->h, l1);
  getHexagonalPixelLineIndex(position, -sqrt(3.), 2.*hp->h, l2);
}



/** Set up the auxiliary array that is used to convert line indices to the corresponding
 * pixel index. The pixel indices start at 0. */
static inline void setLineIndexInformation(HexagonalPixels* hp, int l0, int l1, int l2, 
					   HexagonalPixelLineIndexInformation information)
{
  hp->lineIndexInformation
    [l0+HEXAGONAL_PIXELS_LINE_INDEX_OFFSET]
    [l1+HEXAGONAL_PIXELS_LINE_INDEX_OFFSET]
    [l2+HEXAGONAL_PIXELS_LINE_INDEX_OFFSET] = information;
}



/** Determine the pixel that corresponds to the 3 given line indices. 
 * The function returns information about the sub-triangle, namely the pixel index 
 * and the orientation within the hexagonal pixel it belongs to.
 * If the indices correspond to an invalid pixel, the return value has the pixelindex -1.
 * The pixel indices start at 0.*/
static inline HexagonalPixelLineIndexInformation getLineIndexInformation
(HexagonalPixels* hp, int l0, int l1, int l2)
{
  if ((l0+HEXAGONAL_PIXELS_LINE_INDEX_OFFSET<0) || (l0>HEXAGONAL_PIXELS_LINE_INDEX_OFFSET) ||
      (l1+HEXAGONAL_PIXELS_LINE_INDEX_OFFSET<0) || (l1>HEXAGONAL_PIXELS_LINE_INDEX_OFFSET) ||
      (l2+HEXAGONAL_PIXELS_LINE_INDEX_OFFSET<0) || (l2>HEXAGONAL_PIXELS_LINE_INDEX_OFFSET)) {
    HexagonalPixelLineIndexInformation information = { .pixelindex = INVALID_PIXEL };
    return(information);
  } else {
    return(hp->lineIndexInformation
	   [l0+HEXAGONAL_PIXELS_LINE_INDEX_OFFSET]
	   [l1+HEXAGONAL_PIXELS_LINE_INDEX_OFFSET]
	   [l2+HEXAGONAL_PIXELS_LINE_INDEX_OFFSET]);
  }
}



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
  // Determine the coordinates of the pixel centers (pixel numbering according to Oosterbroek):
  struct Point2d centers[HTRS_N_PIXELS];

  centers[0].x  =  0.; // Pixel number 1
  centers[0].y  =  0.; 
  //
  centers[1].x  =  1.  * hp->h;
  centers[1].y  = -1.5 * hp->a;

  centers[2].x  =  2.  * hp->h;
  centers[2].y  =  0.;

  centers[3].x  =  1.  * hp->h;
  centers[3].y  =  1.5 * hp->a;

  centers[4].x  = -1.  * hp->h;
  centers[4].y  =  1.5 * hp->a;

  centers[5].x  = -2.  * hp->h;
  centers[5].y  =  0.;

  centers[6].x  = -1.  * hp->h;
  centers[6].y  = -1.5 * hp->a;
  //
  centers[7].x  =  0.;
  centers[7].y  = -3.  * hp->a;

  centers[8].x  =  2.  * hp->h;
  centers[8].y  = -3.  * hp->a;

  centers[9].x  =  3.  * hp->h;
  centers[9].y  = -1.5 * hp->a;

  centers[10].x =  4.  * hp->h;
  centers[10].y =  0.;

  centers[11].x =  3.  * hp->h;
  centers[11].y =  1.5 * hp->a;

  centers[12].x =  2.  * hp->h;
  centers[12].y =  3.  * hp->a;

  centers[13].x =  0.;
  centers[13].y =  3.  * hp->a;

  centers[14].x = -2.  * hp->h;
  centers[14].y =  3.  * hp->a;

  centers[15].x = -3.  * hp->h;
  centers[15].y =  1.5 * hp->a;

  centers[16].x = -4.  * hp->h;
  centers[16].y =  0.;

  centers[17].x = -3.  * hp->h;
  centers[17].y = -1.5 * hp->a;

  centers[18].x = -2.  * hp->h;
  centers[18].y = -3.  * hp->a;
  //
  centers[19].x = -1.  * hp->h;
  centers[19].y = -4.5 * hp->a;

  centers[20].x =  1.  * hp->h;
  centers[20].y = -4.5 * hp->a;

  centers[21].x =  3.  * hp->h;
  centers[21].y = -4.5 * hp->a;

  centers[22].x =  4.  * hp->h;
  centers[22].y = -3.  * hp->a;

  centers[23].x =  5.  * hp->h;
  centers[23].y = -1.5 * hp->a;

  centers[24].x =  6.  * hp->h;
  centers[24].y =  0.;

  centers[25].x =  5.  * hp->h;
  centers[25].y =  1.5 * hp->a;

  centers[26].x =  4.  * hp->h;
  centers[26].y =  3.  * hp->a;

  centers[27].x =  3.  * hp->h;
  centers[27].y =  4.5 * hp->a;

  centers[28].x =  1.  * hp->h;
  centers[28].y =  4.5 * hp->a;

  centers[29].x = -1.  * hp->h;
  centers[29].y =  4.5 * hp->a;

  centers[30].x = -3.  * hp->h;
  centers[30].y =  4.5 * hp->a;

  centers[31].x = -4.  * hp->h;
  centers[31].y =  3.  * hp->a;

  centers[32].x = -5.  * hp->h;
  centers[32].y =  1.5 * hp->a;

  centers[33].x = -6.  * hp->h;
  centers[33].y =  0.;

  centers[34].x = -5.  * hp->h;
  centers[34].y = -1.5 * hp->a;

  centers[35].x = -4.  * hp->h;
  centers[35].y = -3.  * hp->a;

  centers[36].x = -3.  * hp->h; // Pixel number 37
  centers[36].y = -4.5 * hp->a;


  // Set up the auxiliary array to convert line indices to pixel indices.
  // First clear the auxiliary array.
  int l0, l1, l2;
  for (l0=0; l0<2*HEXAGONAL_PIXELS_LINE_INDEX_OFFSET+1; l0++) {
    for (l1=0; l1<2*HEXAGONAL_PIXELS_LINE_INDEX_OFFSET+1; l1++) {
      for (l2=0; l2<2*HEXAGONAL_PIXELS_LINE_INDEX_OFFSET+1; l2++) {
	HexagonalPixelLineIndexInformation lineIndexInformation = 
	  { .pixelindex = INVALID_PIXEL };
	hp->lineIndexInformation[l0][l1][l2] = lineIndexInformation;
      }
    }
  }

  // Now determine for each pixel the 6 valid combinations of line indices 
  // corresponding to the 6 individual sub-triangles of this hexagonal pixel.
  int pixel, direction;
  struct Point2d point;
  const double sin60 = sin(M_PI/3.);
  const double cos60 = cos(M_PI/3.);
  for (pixel=0; pixel<HTRS_N_PIXELS; pixel++) {
    // For each pixel choose 6 points located around the center and
    // determine the line indices which define this pixel section.
    for (direction=0; direction<6; direction++) {
      HexagonalPixelLineIndexInformation lineIndexInformation = 
	{ .pixelindex = pixel };

      point = centers[pixel];
	
      switch(direction) {
      case 0: // right from center
	point.x += hp->h/2;
	lineIndexInformation.orientation = TRIANGLE_ORIENTATION_RIGHT;
	break;
      case 1: // upper right section
	point.x += hp->h/2*cos60;
	point.y += hp->h/2*sin60;
	lineIndexInformation.orientation = TRIANGLE_ORIENTATION_UPPER_RIGHT;
	break;
      case 2: // upper left section
	point.x -= hp->h/2*cos60;
	point.y += hp->h/2*sin60;
	lineIndexInformation.orientation = TRIANGLE_ORIENTATION_UPPER_LEFT;
	break;
      case 3: // left from center
	point.x -= hp->h/2;
	lineIndexInformation.orientation = TRIANGLE_ORIENTATION_LEFT;
	break;
      case 4: // lower left section
	point.x -= hp->h/2*cos60;
	point.y -= hp->h/2*sin60;
	lineIndexInformation.orientation = TRIANGLE_ORIENTATION_LOWER_LEFT;
	break;
      case 5: // lower right section
	point.x += hp->h/2*cos60;
	point.y -= hp->h/2*sin60;
	lineIndexInformation.orientation = TRIANGLE_ORIENTATION_LOWER_RIGHT;
	break;
      }

      getHexagonalPixelLineIndices(hp, point, &l0, &l1, &l2);
      setLineIndexInformation(hp, l0, l1, l2, lineIndexInformation);
    } // End of loop over directions
  } // End of loop over all pixels


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



static inline double getHexagonalPixelDistance2Line(struct Point2d position,
						    double m, double t)
{
  return(sqrt( pow(t-position.x+m*position.y, 2.) / (1+pow(m,2.)) ));
}



void getHexagonalPixel(HexagonalPixels* hp, struct Point2d position, int* pixel)
{
  // Determine the 3 lines that define the sub-triangle of the pixel
  // the point position lies within.
  int l0, l1, l2;
  getHexagonalPixelLineIndices(hp, position, &l0, &l1, &l2);

  // From these 3 line indices now determine the pixel.
  HexagonalPixelLineIndexInformation lineIndexInformation = 
    getLineIndexInformation(hp, l0, l1, l2);
  *pixel = lineIndexInformation.pixelindex;
}



int getHexagonalPixelSplits(HexagonalPixels* hp, GenericDetector* gd, 
			    struct Point2d position, 
			    int* pixel, double* fraction)
{
  // Determine the 3 lines that define the sub-triangle of the pixel
  // the point position lies within.
  int l0, l1, l2;
  getHexagonalPixelLineIndices(hp, position, &l0, &l1, &l2);

  // From these 3 line indices now determine the pixel.
  HexagonalPixelLineIndexInformation lineIndexInformation = 
    getLineIndexInformation(hp, l0, l1, l2);
  pixel[0] = lineIndexInformation.pixelindex;


  //    Split Events.
  // Check if charge cloud size is zero.
  if (gd->ccsize < 1.e-20) {
    // Only single events are possible!
    fraction[0] = 1.;
    return(1);
  }


  // Determine the distance to the next neighboring pixel.
  double distance; // Distance to next neighboring pixel
  int k0=l0, k1=l1, k2=l2; // Line indices of next neighboring pixel
  switch (lineIndexInformation.orientation) {
  case TRIANGLE_ORIENTATION_RIGHT:
    distance = getHexagonalPixelDistance2Line(position,        0., (l0+1)  *hp->h);
    k0++;
    break;
  case TRIANGLE_ORIENTATION_UPPER_RIGHT:
    distance = getHexagonalPixelDistance2Line(position, -sqrt(3.), (l2+1)*2*hp->h);
    k2++;
    break;
  case TRIANGLE_ORIENTATION_UPPER_LEFT:
    distance = getHexagonalPixelDistance2Line(position,  sqrt(3.),  l1   *2*hp->h);
    k1--;
    break;
  case TRIANGLE_ORIENTATION_LEFT:
    distance = getHexagonalPixelDistance2Line(position,        0.,  l0     *hp->h);
    k0--;
    break;
  case TRIANGLE_ORIENTATION_LOWER_LEFT:
    distance = getHexagonalPixelDistance2Line(position, -sqrt(3.),  l2   *2*hp->h);
    k2--;
    break;
  case TRIANGLE_ORIENTATION_LOWER_RIGHT:
    distance = getHexagonalPixelDistance2Line(position,  sqrt(3.), (l1+1)*2*hp->h);
    k1++;
    break;
  }

  // Check whether this is really a split event:
  if (distance > gd->ccsize) {
    // Single event!
    fraction[0] = 1.;
    return(1);

  } else {
    // Double event!
    pixel[1] = getLineIndexInformation(hp, k0, k1, k2).pixelindex;

    double mindistgauss = gaussint(distance/gd->ccsigma);
    fraction[0] = 1. - mindistgauss;
    fraction[1] =      mindistgauss;
    
    return(2);
  }
}


