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

#ifndef COMAEVENTFILE_H
#define COMAEVENTFILE_H 1

#include "sixt.h"
#include "comaevent.h"
#include "eventlist.h"


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


typedef struct {

  /** Generic EventFile data structure. */
  EventList generic;

  /* Column numbers of the individual Coded mask specific event list
      entries. The numbers start at 1. The number 0 means, that there
      is no corresponding column in the table. */
  int ctime, ccharge, crawx, crawy;

} CoMaEventFile;


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////


/** Opens an existing FITS file with a binary table event list. Apart
    from opening the FITS file the function also determines the number
    of rows in the FITS table and initializes the CoMaEventFile data
    structure. The access_mode parameter can be either READONLY or
    READWRITE. */
CoMaEventFile* openCoMaEventFile(char* const filename,
				 const int access_mode,
				 int* const status);

/** Create and open a new FITS event file for the HTRS detector from a
    given FITS template. If the file already exists, the old file is
    deleted and replaced by an empty one. Apart from opening the FITS
    file the function also initializes the CoMaEventFile data
    structure by calling openCoMaEventFile(). The access_mode
    parameter is always READWRITE. */
CoMaEventFile* openNewCoMaEventFile(char* const filename,
				    char* const template,
				    int* const status);

/** Close an open Coded Mask event list FITS file. */
int closeCoMaEventFile(CoMaEventFile* ef);

/** Append a new CoMaEvent to the to event list. In the given
    CoMaEvent data structure the pixel numbering starts at 0, but in
    the event file the numbering has to start at 1. So the routine
    adds a 1 to the pixel index. The return value is the error
    status. */
int addCoMaEvent2File(CoMaEventFile* ef, CoMaEvent* event);

/** Read the next CoMaEvent from the CoMaEventFile. This routine
    increases the internal counter of the CoMaEventFile data
    structure. In the event file the numbering of the pixels starts
    at 1, whereas in the returned CoMaEvent data structure the
    numbering starts at 0. The return value is the error status. */
int CoMaEventFile_getNextRow(CoMaEventFile* ef, CoMaEvent* event);


#endif /* COMAEVENTFILE */
