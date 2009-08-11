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
static inline void setLineIndices2Pixels(HexagonalPixels* hp, int l0, int l1, int l2, int pixel)
{
  hp->LineIndices2Pixel
    [l0+HEXAGONAL_PIXELS_LINE_INDEX_OFFSET]
    [l1+HEXAGONAL_PIXELS_LINE_INDEX_OFFSET]
    [l2+HEXAGONAL_PIXELS_LINE_INDEX_OFFSET] = pixel;
}



/** Determine the pixel that corresponds to the 3 given line indices. 
 * The pixel indices start at 0. */
static inline int getHexagonalPixelFromLineIndices(HexagonalPixels* hp, int l0, int l1, int l2)
{
  if ((l0+HEXAGONAL_PIXELS_LINE_INDEX_OFFSET<0) || (l0>HEXAGONAL_PIXELS_LINE_INDEX_OFFSET) ||
      (l1+HEXAGONAL_PIXELS_LINE_INDEX_OFFSET<0) || (l1>HEXAGONAL_PIXELS_LINE_INDEX_OFFSET) ||
      (l2+HEXAGONAL_PIXELS_LINE_INDEX_OFFSET<0) || (l2>HEXAGONAL_PIXELS_LINE_INDEX_OFFSET)) {
    return(INVALID_PIXEL);
  } else {
    return(hp->LineIndices2Pixel
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
	hp->LineIndices2Pixel[l0][l1][l2] = INVALID_PIXEL;
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
      point = centers[pixel];
	
      switch(direction) {
      case 0: // right from center
	point.x += hp->h/2;
	break;
      case 1: // upper right section
	point.x += hp->h/2*cos60;
	point.y += hp->h/2*sin60;
	break;
      case 2: // upper left section
	point.x -= hp->h/2*cos60;
	point.y += hp->h/2*sin60;
	break;
      case 3: // left from center
	point.x -= hp->h/2;
	break;
      case 4: // lower left section
	point.x -= hp->h/2*cos60;
	point.y -= hp->h/2*sin60;
	break;
      case 5: // lower right section
	point.x += hp->h/2*cos60;
	point.y -= hp->h/2*sin60;
	break;
      }

      getHexagonalPixelLineIndices(hp, point, &l0, &l1, &l2);
      setLineIndices2Pixels(hp, l0, l1, l2, pixel);
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



/** Determines the minimum distance value out of an array with 6
 * entries and returns the corresponding index. */
static inline int getMinimumDistance(double array[]) 
{
  int count, index=0;
  double minimum=array[0];

  for (count=1; count<6; count++) {
    if ( (minimum < 0.) ||
	 ((array[count]<=minimum)&&(array[count]>=0.)) ) {
      minimum = array[count];
      index = count;
    }
  }

  return(index);
}



///////////////////////////////////////////
static inline double getHexagonalPixelDistance2Line(struct Point2d position,
						    double m, double t)
{
  return(sqrt( pow(t-position.x+m*position.y, 2.) / (1+pow(m,2.)) ));
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
  pixel[0] = getHexagonalPixelFromLineIndices(hp, l0, l1, l2);


  //    Split Events.
  // Check if charge cloud size is zero.
  if (gd->ccsize < 1.e-20) {
    // Only single events are possible!
    fraction[0] = 1.;
    return(1);
  }

  // Distances to neighbouring pixel segments (equilateral triangles, CAN possibly
  // belong to the same pixel).
  double distances[6] = {
    // right
    getHexagonalPixelDistance2Line(position,        0., (l0+1)  *hp->h),
    // upper right
    getHexagonalPixelDistance2Line(position, -sqrt(3.), (l2+1)*2*hp->h),
    // upper left
    getHexagonalPixelDistance2Line(position,  sqrt(3.),  l1   *2*hp->h),
    // left
    getHexagonalPixelDistance2Line(position,        0.,  l0     *hp->h),
    // lower left
    getHexagonalPixelDistance2Line(position, -sqrt(3.),  l2   *2*hp->h),
    // lower right
    getHexagonalPixelDistance2Line(position,  sqrt(3.), (l1+1)*2*hp->h),
  };

  // Find the closest distance to the nearest neighbouring pixel.
  int count, mindist;
  double minimum;
  int dl[3][6] = {
    {1, 0, 0, -1, 0, 0},
    {0, 0, -1, 0, 0, 1},
    {0, 1, 0, 0, -1, 0}
  };
  for(count=0; count<3; count++) {
    mindist = getMinimumDistance(distances);
    minimum = distances[mindist];

    int k0 = l0+dl[0][mindist];
    int k1 = l1+dl[1][mindist];
    int k2 = l2+dl[2][mindist];
    pixel[1] = getHexagonalPixelFromLineIndices(hp, k0, k1, k2);
    
    if (pixel[1] != pixel[0]) break;
    distances[mindist] = -1.;
  }

  if ((minimum > gd->ccsize) || (pixel[1] == pixel[0])) {
    // Single event!
    fraction[0] = 1.;
    return(1);

  } else {
    // Double event!
    double mindistgauss = gaussint(minimum/gd->ccsigma);
    fraction[0] = 1. - mindistgauss;
    fraction[1] =      mindistgauss;
    
    return(2);
  }
}


