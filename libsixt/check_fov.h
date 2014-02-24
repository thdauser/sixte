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


   Copyright 2007-2014 Christian Schmid, FAU
*/

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
	      const double min_align);


#endif /* CHECK_FOV_H */
