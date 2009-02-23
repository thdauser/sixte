#ifndef CHECK_FOV_H
#define CHECK_FOV_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vector.h"


int check_fov(const struct vector x, const struct vector x0, /* struct vector h1, struct vector h2, double sin_dec_max, double sin_rasc_max, */ const double min_dev);

#endif
