#include "check_fov.h"

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
inline int check_fov(Vector* const x, Vector* const x0, 
		     /* Vector h1, Vector h2, double sin_dec_max, 
			double sin_rasc_max, */ 
		     const double min_align)
{
  int result = 0; // source is inside the FOV !

  if (scalar_product(x, x0) < min_align) {
    result = 3;   // source is outside the FOV !
  }

  /* else {
    if (fabs(scalar_product(x,h1)) > sin_dec_max) {
      result += 1;
    }
    if (fabs(scalar_product(x,h2)) > sin_rasc_max) {
      result += 2;
    }
    /////////////////
    // probably it is not necessary to check right ascension, 
    // if declination is already out of range !!
    /////////////////
  } */

  return(result);
}


