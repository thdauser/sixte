/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

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
