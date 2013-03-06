#include "det_phi_max.h"

float det_phi_max(const float d,
		  const float x_mask,
		  const float y_mask,
		  const float x_det,
		  const float y_det)
{
  float phi_max;
  phi_max = atan(1/(2*d)*sqrt((x_mask+x_det)*(x_mask+x_det)+(y_mask+y_det)*(y_mask+y_det)));

  return(phi_max);
}
