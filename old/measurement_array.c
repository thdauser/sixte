#include "measurement_array.h"


/////////////////////////////////////////////////////////////////
// Function plots the content of a square array to a png file.
int plot_array(double ** array, int width, int plotpixelwidth, char filename[]) {
  int x, y;             // counters
  double maximum = 0.;  // value of the brightest point in the array

  // variables for png-access:
  bitmap_t bitmap;
  pixel_t pixel;

  // initialize image:
  bitmap.width = width*plotpixelwidth;
  bitmap.height = width*plotpixelwidth;
  bitmap.name = filename;

  if(bmp_init(&bitmap) !=IMGLIB_OK){
    printf("Bitmap initialization failed\n");
    return(IMGLIB_ERROR);
  }


  // find the brightest pixel in the array (to scale the brightness of the picture):
  for(x=0; x < width; x++) {
    for(y=0; y < width; y++) {
      if(array[x][y] > maximum) {
        maximum = array[x][y];
      }
    }
  }


  // create the png-plot:
  for(x=0; x < bitmap.height; x++) {
    for(y=0; y < bitmap.width; y++) {
      pixel.r = (unsigned char)
	(255.0*array[x/plotpixelwidth][y/plotpixelwidth]/maximum);
      pixel.g = 0;
      pixel.b = 0;

      if(bmp_set_pixel_rgb(&bitmap, &pixel, x, y)){
        printf("Failed to set RGB of pixel (%d,%d)\n",x,y);
        return(IMGLIB_ERROR);
      }
    }
  }

  if(png_write(&bitmap)!=IMGLIB_OK){
    printf("PNG output failed\n");
    return(IMGLIB_ERROR);
  }

  return(IMGLIB_OK);
}



///////////////////////////////////////////////
// Function clears the given array.
void clear_array(float **array, int width) {
  int x,y;

  for (x=0; x<width; x++) {
    for (y=0; y<width; y++) {
      array[x][y] = 0.;
    }
  }
}


