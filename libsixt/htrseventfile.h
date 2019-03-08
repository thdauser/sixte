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

#ifndef HTRSEVENTFILE_H
#define HTRSEVENTFILE_H 1

#include "sixt.h"
#include "htrsevent.h"
#include "eventlist.h"


typedef struct {

  /** Generic EventList data structure. */
  EventList generic;

  /* Column numbers of the individual HTRS-specific event list entries.
   * The numbers start at 1. The number 0 means, that there
   * is no corresponding column in the table. */
  int ctime, cpha, cenergy, cpixel, cgrade, cx, cy;

} HTRSEventFile;


/////////////////////////////////////////////////////////////////////


/** Opens an existing FITS file with a binary table event list. Apart
    from opening the FITS file the function also determines the number
    of rows in the FITS table and initializes the HTRSEventFile data
    structure. The access_mode parameter can be either READONLY or
    READWRITE. */
int openHTRSEventFile(HTRSEventFile* hef,
		      char* const filename,
		      const int access_mode);

/** Create and open a new FITS event file for the HTRS detector from a
    given FITS template. If the file already exists, the old file is
    deleted and replaced by an empty one.  Apart from opening the FITS
    file the function also initializes the HTRSEventFile data
    structure by calling openHTRSEventFile(). The access_mode
    parameter is always READWRITE. */
int openNewHTRSEventFile(HTRSEventFile* hef,
			 char* const filename,
			 char* const template);

/** Close an open HTRS event list FITS file. */
int closeHTRSEventFile(HTRSEventFile* hef);

/** Append a new HTRS event to the to event list. In the given
    HTRSEvent data structure the pixel numbering starts at 0, but in
    the event file the numbering has to start at 1. So the routine
    adds a 1 to the pixel index. The return value is the error
    status. */
int addHTRSEvent2File(HTRSEventFile* hef, HTRSEvent* event);

/** Read a specific row from the HTRSEventfile. This routine does NOT
    increase the internal counter of the HTRSEventFile data structure.
    In the event file the numbering of the pixels starts at 1, whereas
    in the returned HTRSEvent data structure the numbering starts at
    0. The return value of the function is the error status. */
int HTRSEventFile_getRow(HTRSEventFile* hef, HTRSEvent* event, const long row);

/** Read the next HTRSEvent from the HTRSEventFile. This routine
    increases the internal counter of the HTRSEventFile data
    structure. In the event file the numbering of the pixels starts at
    1, whereas in the returned HTRSEvent data structure the numbering
    starts at 0. The return value is the error status. The return
    value is the error status. */
int HTRSEventFile_getNextRow(HTRSEventFile* hef, HTRSEvent* event);

/** Write the specified column in the HTRSEventFile. The data to be
    written is given by the HTRSEvent. The row already must
    exist. Otherwise an error is returned. */
int HTRSEventFile_writeRow(HTRSEventFile*, HTRSEvent*, const long row);


#endif /* HTRSEVENTFILE */
