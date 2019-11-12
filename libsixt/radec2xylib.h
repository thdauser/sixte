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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef RADEC2XYLIB_H
#define RADEC2XYLIB_H 1

#include "event.h"
#include "wcs.h"
#include "eventfile.h"

////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////

/** Struct holding image x, z positions */
typedef struct {
    long x, y;
} ImgPos;


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** */
ImgPos radec2xy( Event* const event, struct wcsprm* const wcs, int* const status );

struct wcsprm getRadec2xyWCS( float* RefRA, float* RefDec, char* Projection, int* const status );

void addXY2eventfile(EventFile* const evtfile, float* RefRA, float* RefDec, char* Projection, int* const status);

#endif /* RADEC2XY_H */
