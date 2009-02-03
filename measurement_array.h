#ifndef MEASUREMENT_ARRAY_H
#define MEASUREMENT_ARRAY_H 1

#include <png.h>

#include "imglib.h"


// clears the given array:
void clear_array(float **array, int width);

// plots the content of a square array to a png file.
int plot_array(double **array, int width, int plotpixelwidth, char filename[]);


#endif
