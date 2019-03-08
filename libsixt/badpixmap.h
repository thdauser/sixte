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

#ifndef BADPIXMAP_H
#define BADPIXMAP_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Constants.
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Bad pixel map. */
typedef struct {

  /** Dimensions of the bad pixel map. */
  int xwidth, ywidth;

  /** 2-dimensional pixel array containing the bad pixel map. */
  float** pixels;

  /** This flag specifies if the column contains any bad pixels (value
      1). If not (value 0), the bad pixel application routine can jump
      to the next column. */
  int* anybadpix;

} BadPixMap;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Allocates memory for a new empty BadPixMap data
    structure. */
BadPixMap* newBadPixMap(int* const status);

/** Allocate memory for a new BadPixMap data structure via the
    constructor and load the bad pixel data from the specified
    file. */
BadPixMap* loadBadPixMap(const char* const filename, int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the BadPixMap data structure to NULL. */
void destroyBadPixMap(BadPixMap** const map);


#endif /* BADPIXMAP_H */
