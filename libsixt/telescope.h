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

#ifndef TELESCOPE_H
#define TELESCOPE_H 1

#include "vector.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


struct Telescope {
  /* Unit vectors defining the telescope attitude. */

  /** Direction of the telescope axis (in simplest case same direction
      as r). */
  Vector nz;
  /** Usually points in the direction of motion, but is perpendicular
      to nz. */
  Vector nx;
  /** Perpendicular to nx and nz */
  Vector ny; 
};


#endif /* TELESCOPE_H */

