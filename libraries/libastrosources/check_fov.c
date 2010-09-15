#include "check_fov.h"


int check_fov(const Vector* const x, const Vector* const x0, 
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


