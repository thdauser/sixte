/* imglib.c - A library to initialise a bitmap, set each pixel's RGB and write
 *            out to a PNG compressed file
 *
 * Copyright (C) 2005 - 2006
 *                    Joern Wilms
 *                    Robert Dowse
 *                    Robert Swain
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <png.h>

#include "imglib.h"
#include "print_errno_info.c"

int bmp_init(bitmap_t *bitmap) {

  bitmap->bmp=NULL;
  if((bitmap->height > 0) && (bitmap->width > 0)) {
    bitmap->bmp=((unsigned char *)malloc(bitmap->width*bitmap->height*3*sizeof(unsigned char)));
  }

  if (bitmap->bmp!=NULL) return(IMGLIB_OK);
  return(IMGLIB_ERROR);

}

int bmp_set_pixel_rgb(bitmap_t *bitmap, pixel_t *pixel, int bmp_x, int bmp_y) {
  if ( (bmp_x<0) || (bmp_x>=bitmap->width) ||
       (bmp_y<0) || (bmp_y>=bitmap->height)) {
     return(IMGLIB_ERROR);
  }

  // byte offset of pixel from bmp[]
  // i.e. y rows of width rgb pixels + x rgb pixels
  int pixel_pos = (bmp_y*bitmap->width + bmp_x)*3;

  // assign r g b values accordingly
  bitmap->bmp[pixel_pos] = pixel->r;
  bitmap->bmp[pixel_pos + 1] = pixel->g;
  bitmap->bmp[pixel_pos + 2] = pixel->b;
  return(IMGLIB_OK);

}

int bmp_set_pixel_grey(bitmap_t *bitmap, unsigned char grey, int bmp_x, int bmp_y) {
  pixel_t greyval;
  greyval.r=grey;
  greyval.g=grey;
  greyval.b=grey;

  return bmp_set_pixel_rgb(bitmap,&greyval,bmp_x,bmp_y);
}

int png_write(bitmap_t *bitmap) {

// PNG writing using libpng and the array of pointers to rows method

  int i;
  unsigned char *row_ptrs[bitmap->height];

  FILE *fp = fopen(bitmap->name, "wb");

  if (!fp) {
    print_errno_info(errno);
    return (IMGLIB_ERROR); // error if the file did not open correctly
  }

  // calulate the pointers to each row of the image for use later
  for(i = 0; i < bitmap->height; i++) {
    row_ptrs[bitmap->height-i-1] = bitmap->bmp + bitmap->width*3*i;
  }

  // allocate and initialise the png struct
  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr)
    return (IMGLIB_ERROR);

  // allocate and initialise the info struct
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr)
  {
    png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
    return (IMGLIB_ERROR);
  }

  // initialise the io
  if (setjmp(png_jmpbuf(png_ptr))) {
    fprintf(stderr,"Error during png_init_io\n");
    return(IMGLIB_ERROR);
  }
  png_init_io(png_ptr, fp);


  // set the png header info
  if (setjmp(png_jmpbuf(png_ptr))) {
    fprintf(stderr,"Error during png_set_IHDR\n");
    return(IMGLIB_ERROR);
  }
  png_set_IHDR(
    png_ptr, // pointer to the png file
    info_ptr, // pointer to the header info within said png file
    bitmap->width,
    bitmap->height,
    8, // bit depth per channel
    PNG_COLOR_TYPE_RGB, // rgb colour type
    PNG_INTERLACE_NONE, // no interlacing
    PNG_COMPRESSION_TYPE_DEFAULT,
    PNG_FILTER_TYPE_DEFAULT);

  // write the header info
  if (setjmp(png_jmpbuf(png_ptr))) {
    fprintf(stderr,"Error during png_write_info\n");
    return(IMGLIB_ERROR);
  }
  png_write_info(png_ptr, info_ptr);
  // write the image using the row pointers
  if (setjmp(png_jmpbuf(png_ptr))) {
    fprintf(stderr,"Error during png_write_image\n");
    return(IMGLIB_ERROR);
  }
  png_write_image(png_ptr, row_ptrs);
  // finish writing and free allocated memory
  if (setjmp(png_jmpbuf(png_ptr))) {
    fprintf(stderr,"Error during png_write_end\n");
    return(IMGLIB_ERROR);
  }
  png_write_end(png_ptr, info_ptr);

  if (setjmp(png_jmpbuf(png_ptr))) {
    fprintf(stderr,"Error during png_destroy_write_struct\n");
    return(IMGLIB_ERROR);
  }
  png_destroy_write_struct(&png_ptr, &info_ptr);


  // close picture output-file
  fclose(fp);

  return(IMGLIB_OK);

}
