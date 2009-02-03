#include <stdio.h>
#include <stdlib.h>
#include <math.h>



#define INVALID_PIXEL -1

#define DET_WIDTH 7      // maximum number of pixels in one line
#define N_PIXELS  37     // total number of pixels


struct Point2i {
  int x, y;
};

struct Point2d {
  double x, y;
};



const double a     = 1.0;          // length of a pixel edge
const double sqrt3 = 1.732050808;  // sqrt(3.);
const double h = 0.5 * 1.0 * 1.732050808;
//                     |     |-> sqrt(3.)
//                     |-> a



///////////////////////////////////////////
double linear_function(double x, double m, double t)
{
  return(m*x + t);
}



///////////////////////////////////////////
// Returns the lower line index of one particular line group with specified m
// and variable t for a given point in the 2D plane.
int get_line(struct Point2d point, double m, double dt)
{
  if (point.y < linear_function(point.x, m, -DET_WIDTH*dt)) {
    return(INVALID_PIXEL);
  } else if (point.y > linear_function(point.x, m, DET_WIDTH*dt)) {
    return(INVALID_PIXEL);
  } else {
    int counter;
    for (counter=0; counter<2*DET_WIDTH; counter++) {
      if (point.y < linear_function(point.x, m, (counter+1-DET_WIDTH)*dt))
	break;
    }
    return(counter);
  }
}




///////////////////////////////////////////
void get_lines(struct Point2d point, int* l0, int* l1, int* l2)
{
  *l0 = get_line(point,     0.,   h);
  *l1 = get_line(point, -sqrt3, 2*h);
  *l2 = get_line(point,  sqrt3, 2*h);
}




///////////////////////////////////////////
// This function determines the pixel index for a given 2D floating point.
int get_pixel(int*** pixel_relations, struct Point2d point)
{
  int l0, l1, l2;
  get_lines(point, &l0, &l1, &l2);

  if ((l0==INVALID_PIXEL)||(l1==INVALID_PIXEL)||(l2==INVALID_PIXEL)) {
    return(INVALID_PIXEL);
  } else {
    return(pixel_relations[l0][l1][l2]);
  }
}




///////////////////////////////////////////
void initialize_pixel_relations(int**** pixel_relations, struct Point2d* centers)
{
  int l0, l1, l2;

  // Get memory and clear the array
  *pixel_relations = (int***) malloc(2*DET_WIDTH * sizeof(int**));
  for (l0=0; l0<2*DET_WIDTH; l0++) {
    (*pixel_relations)[l0] = (int**) malloc(2*DET_WIDTH * sizeof(int*));
    for (l1=0; l1<2*DET_WIDTH; l1++) {
      (*pixel_relations)[l0][l1] = (int*) malloc(2*DET_WIDTH * sizeof(int));
      for (l2=0; l2<2*DET_WIDTH; l2++) {
	(*pixel_relations)[l0][l1][l2] = INVALID_PIXEL;
      }
    }
  }

  int pixel, direction;
  struct Point2d point;
  const double sin30 = sin(M_PI/6.);
  const double cos30 = cos(M_PI/6.);

  for(pixel=0; pixel<N_PIXELS; pixel++) {
    // For each pixel choose 6 points located around the center and
    // determine the line indices which define this pixel section.
    
    for(direction=0; direction<6; direction++) {
      point    = centers[pixel];

      switch(direction) {
      case 0: // above center
	point.y += h/2;
	break;
      case 1: // upper right section
	point.x += h/2*cos30;
	point.y += h/2*sin30;
	break;
      case 2: // lower right section
	point.x += h/2*cos30;
	point.y -= h/2*sin30;
	break;
      case 3: // below center
	point.y -= h/2;
	break;
      case 4: // lower left section
	point.x -= h/2*cos30;
	point.y -= h/2*sin30;
	break;
      case 5: // upper left section
	point.x -= h/2*cos30;
	point.y += h/2*sin30;
	break;
      }

      get_lines(point, &l0, &l1, &l2);
      (*pixel_relations)[l0][l1][l2] = pixel;
    }
  }

}




///////////////////////////////////////////
int main()
{
  // -- 2 different access numbering schemes --
  // 2d array contains linear numbers
  int icoordinates2pixel[DET_WIDTH][DET_WIDTH];   
  // linear array contains 2d integer coordinates
  struct Point2i pixel2icoordinates[N_PIXELS]; 

  struct Point2d centers[N_PIXELS];      // centers of the HTRS pixels

  // Data structure to obtain a pixel from a given coordinate:
  int*** lines2pixel; 
 
  // Counters
  int xi, yi;
  int pixel_counter;



  // Set the relation between the two different numbering arrays of 
  // the pixels in the hexagonal structure.
  pixel_counter=0;
  for(xi=0; xi<DET_WIDTH; xi++) {
    for(yi=0 ; yi<DET_WIDTH; yi++) {
      if (((xi==0)||(xi==6)) && ((yi==0)||(yi>=5))) {
	icoordinates2pixel[xi][yi] = INVALID_PIXEL;
      } else if (((xi==1)||(xi==5)) && ((yi==0)||(yi==6))) {
	icoordinates2pixel[xi][yi] = INVALID_PIXEL;
      } else if (((xi==2)||(xi==4)) && (yi==6)) {
	icoordinates2pixel[xi][yi] = INVALID_PIXEL;
      } else {
	icoordinates2pixel[xi][yi] = pixel_counter;
	pixel2icoordinates[pixel_counter].x = xi;
	pixel2icoordinates[pixel_counter].y = yi;
	pixel_counter++;
      }
    }
  }
  // Now the 2 different numbering schemes can be easily converted among each other.


  // Calculate the centers of the hexagonal HTRS pixels.
  for(xi=0; xi<DET_WIDTH; xi++) {
    for(yi=0 ; yi<DET_WIDTH; yi++) {
      pixel_counter = icoordinates2pixel[xi][yi];
      if (pixel_counter != INVALID_PIXEL) {
	centers[pixel_counter].x = ( xi-3 )                  *1.5*a;
	centers[pixel_counter].y = ( yi-3 + ((xi+1)%2)*0.5 ) *2.0*h;
      }
    }
  }


  // Initialize the pixel relations array
  initialize_pixel_relations(&lines2pixel, centers);


  struct Point2d test;
  test.x = 1.5*a;
  test.y = 5.9*h;
  int pixel = get_pixel(lines2pixel, test);
  printf("pixel %d\n", pixel);
  if (pixel != INVALID_PIXEL) 
    printf("coordinates %d %d\n", pixel2icoordinates[pixel].x, 
	   pixel2icoordinates[pixel].y);

  /*
  // Print the pixel centers
  for(pixel_counter=0; pixel_counter<N_PIXELS; pixel_counter++) {
    printf("%lf %lf\n", centers[pixel_counter].x, centers[pixel_counter].y);
  }
  */


  return(EXIT_SUCCESS);
}

