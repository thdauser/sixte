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

#ifndef ATTITUDEFILE_H
#define ATTITUDEFILE_H 1

#include "sixt.h"


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** Data structure describing an attitude file. */
typedef struct {
  /** File pointer to the FITS file. */
  fitsfile* fptr; 

  /** Current row in the attitude FITS file (starting at 0). */
  long row;

  /** Total number of rows in the attitude file. */
  long nrows; 

  /* Column numbers of the individual attitude file entries. The
     numbers start at 1. The number 0 means, that there is no
     corresponding column in the table. */
  int ctime, cra, cdec, crollang;

} AttitudeFile;


/** Contains a line of attitude data from a FITS file. */
typedef struct {
  /** Time for which the AttitudeFileEntry is valid. */
  double time; 
  /** Right ascension of telescope pointing direction [deg]. */
  float ra; 
  /** Declination of telescope pointing direction [deg].*/
  float dec; 
  /** Rollangle [deg]. */
  float rollang; 

} AttitudeFileEntry;


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////


/** Reads a line of data from the attitude table in a FITS file. The
    routine does NOT increment the row counter of the AttitudeFile
    object. */
AttitudeFileEntry read_AttitudeFileEntry(AttitudeFile* const af, int* const status);

/** Opens an existing attitude file. The access_mode parameter can be
    either READONLY or READWRITE. */
AttitudeFile* open_AttitudeFile(const char filename[], const int access_mode, 
				int* const status);


#endif /* ATTITUDEFILE_H */
