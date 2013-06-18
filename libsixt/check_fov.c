#include "check_fov.h"


int check_fov(const Vector* const x, const Vector* const x0, 
	      const double min_align)
{
  // By default assume that the source is inside the FOV.
  int result = 0; 

  if (scalar_product(x, x0) < min_align) {
    // The source is located outside the FOV.
    result = 3;
  }

  return(result);
}


