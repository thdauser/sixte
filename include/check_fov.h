#ifndef CHECK_FOV_H
#define CHECK_FOV_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vector.h"


inline int check_fov(Vector* const x, Vector* const x0, 
		     /* Vector h1, Vector h2, double sin_dec_max, 
			double sin_rasc_max, */ 
		     const double min_align);

#endif /* CHECK_FOV_H */
