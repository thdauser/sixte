#ifndef CHECK_FOV_H
#define CHECK_FOV_H 1

#include "sixt.h"
#include "vector.h"


/** Checks, whether the specified source x is inside the FOV.  The
    first check is, whether the source is inside a circle, by 

    <x,x0> < MIN_VAL.  

    The second and more closer check is carried out by the following
    formula:
 
  sin(theta_max) >= fabs(sin(theta)) = fabs(<r,h1>) = fabs(<x-x0,h1>) = 
  fabs(<x,h1> - diff1) = fabs(<x,h1>)
  (diff<i> = 0)
 
    Analogously with phi, h2 and diff2.
 
    The return value is 0, if the source is inside the FOV and >0 if
    it is outside. */
int check_fov(const Vector* const x, const Vector* const x0, 
	      /* Vector h1, Vector h2, double sin_dec_max, 
		 double sin_rasc_max, */ 
	      const double min_align);

#endif /* CHECK_FOV_H */
