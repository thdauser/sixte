/* imglib.h - Useful structs and definitions for imglib
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

#ifndef _IMGLIB_H
#define _IMGLIB_H

// the image
typedef struct {
  int width, height;
  char *name;
  unsigned char *bmp;
} bitmap_t;

// color value for one pixel
typedef struct {
  unsigned char r, g, b;
} pixel_t;

int bmp_init(bitmap_t *bitmap);
int bmp_set_pixel_rgb(bitmap_t *bitmap,pixel_t *pixel, 
		      int bmp_x, int bmp_y);

int bmp_set_pixel_grey(bitmap_t *bitmap,unsigned char grey, 
		       int bmp_x, int bmp_y);

int png_write(bitmap_t *bitmap);

#define IMGLIB_OK 0
#define IMGLIB_ERROR 1


#endif
