#ifndef CHECK_FOV_H
#define CHECK_FOV_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vector.h"


inline int check_fov(struct vector* const x, struct vector* const x0, 
		     /* struct vector h1, struct vector h2, double sin_dec_max, 
			double sin_rasc_max, */ 
		     const double min_align);

#endif /* CHECK_FOV_H */
